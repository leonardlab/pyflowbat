# [ ] inputs for functions into **kwargs
import os
import FlowCal as fc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statsmodels.api as sm

from . import _std_vals

class Workspace:

    def __init__(self, stylesheet = _std_vals.std_pfb_style, lims_file = "_std", full_output = False) -> None:
        self.full_output = full_output
        if lims_file == "_std":
            self.lims = _std_vals.std_lims
        else:
            self.lims = self._read_lims_file(lims_file)
        self.conversion_factors = None
        self.flow_data = None
        self.flow_statistics = None
        self.sample_collections = {}
        self.stats_collections = {}
        self.compensation_matrix = None
        mpl.rcParams["figure.dpi"] = 150
        if stylesheet is not None:
            mpl.style.use(stylesheet)
        #[ ] add R functionality
        self.r_ready = False

    def _verify_R_installation(self):
        from subprocess import CalledProcessError
        from subprocess import check_call
        try:
            check_call(['which', 'R'])
        except CalledProcessError:
            raise RuntimeError("No R installation could be found")
        return True

    def init_r(self, check_R_installation = True):
        if self.r_ready:
            print("R functionality already initialized")
            return
        if check_R_installation:
            self._verify_R_installation()
        from . import r_gating
        self.r_ready = True

    def _read_lims_file(self, lims_file) -> dict[str, list[int]]:
        lims_data = pd.read_csv(lims_file)
        lims_dict = {}
        for header in list(lims_data.columns):
            lims_dict[header] = list(lims_data[header])
        return lims_dict

    def _read_bead_conversion_file(self, conversions_file) -> dict[str, np.ndarray]:
        beads_data = pd.read_csv(conversions_file)
        beads_dict = {}
        for header in list(beads_data.columns):
            beads_dict[header] = np.asarray(beads_data[header])
        return beads_dict

    def _perform_beads_calculations(self, beads_file, beads_fluorescent_channels, beads_num_pops, conversions_file) -> dict[str, list[float, float]]:
        if conversions_file == "_std":
            conversions_data = _std_vals.std_beads_conversions
        else:
            conversions_data = self._read_bead_conversion_file(conversions_file)
        beads_data = fc.io.FCSData(beads_file)
        pops = fc.mef.clustering_gmm(
            beads_data[:, [beads_fluorescent_channels[0][0]]],
            beads_num_pops,
            tol=1e-07,
            min_covar=None,
            scale="logicle",
        )
        pops = np.asarray(pops)
        num_chs = len(beads_fluorescent_channels)
        beads_means = np.zeros((num_chs, beads_num_pops))
        for pop in range(beads_num_pops):
            for count, ch in enumerate(beads_fluorescent_channels):
                tmp_mean = np.mean((beads_data[:, ch[0]])[np.where(pops == pop)])
                beads_means[count][pop] = tmp_mean
        models = [0] * len(beads_fluorescent_channels)
        for count, ch in enumerate(beads_fluorescent_channels):
            models[count] = sm.OLS(endog=conversions_data[ch[1]], exog=beads_means[count]).fit()
            if self.full_output:
                print(models[count].summary())
                plt.figure()
                plt.scatter(beads_means[count], conversions_data[ch[1]])
                plt.plot(
                    [0, np.max(beads_means[count])],
                    [0, np.max(beads_means[count]) * models[count].params[0]],
                )
                plt.show()
        conv_fctrs = {}
        for count, ch in enumerate(beads_fluorescent_channels):
            conv_fctrs[ch[0]] = models[count].params[0]
            conv_fctrs["" + ch[0] +"_stderr"] = models[count].bse[0]
        return conv_fctrs

    def calculate_beads_factors(self, beads_file, beads_fluorescent_channels, beads_num_pops, beads_conversions_file = "_std"):
        self.conversion_factors = self._perform_beads_calculations(beads_file, beads_fluorescent_channels, beads_num_pops, beads_conversions_file)

    def _load_samples(self, file_folder, file_quals):
        extracted_data = {}
        file_names = os.listdir(file_folder)
        file_names.sort()
        for i in range(len(file_names)):
            file_name = str(file_names[i])
            if file_quals(file_name):
                file_path = os.path.join(file_folder, file_names[i])
                fcs_data = fc.io.FCSData(file_path)
                extracted_data[file_names[i]] = fcs_data
        return extracted_data
    
    def load_samples(self, sample_collection_name, samples_folder, samples_quals):
        self.sample_collections[sample_collection_name] = self._load_samples(samples_folder, samples_quals)

    def _extract_statistics(self, data, reqs, columns):
        columns_list = [""] * len(columns)
        for i in range(len(columns)):
            columns_list[i] = columns[i][0]
        destination = pd.DataFrame(columns=columns_list)
        data_names = list(data.copy().keys())
        for i in range(len(data_names)):
            file_name = str(data_names[i])
            if reqs(file_name):
                fcs_data = data[file_name]
                row = [None] * len(columns)
                for j in range(len(columns)):
                    row[j] = columns[j][1](file_name, fcs_data)
                destination.loc[i] = row
        return destination

    def extract_statistics(self, sample_collection_name, statistics_collection_name, samples_quals, statistics_columns):
        self.stats_collections[statistics_collection_name] = self._extract_statistics(self.sample_collections[sample_collection_name], samples_quals, statistics_columns)

    def combine_replicates(self, statistics_collection_name, combined_statistics_collection_name, replicate_definition, columns):
        #[ ]: vectorize
        data = self.stats_collections[statistics_collection_name].copy()
        num_errs = 1
        for i in range(len(columns)):
            if columns[i][2] == True:
                num_errs += 1
        num_columns = len(columns) + num_errs
        columns_list = [""] * num_columns
        num_errs = 0
        for i in range(len(columns)):
            columns_list[i] = columns[i][0]
            if columns[i][2] == True:
                columns_list[len(columns) + num_errs] = columns[i][0] + "_stdErr"
                num_errs += 1
        destination = pd.DataFrame(columns=columns_list)
        k = 0
        i = 0
        j = 0
        while i < (len(data) - 1):
            start_repl = replicate_definition(data.iloc[i])
            new_data = [""] * num_columns
            num_errs = 0
            for l in range(len(columns)):
                j = i + 1
                if columns[l][2] == True:
                    new_data[l] = [columns[l][1](data.iloc[i])]
                    while j < len(data) and start_repl == replicate_definition(data.iloc[j]):
                        new_data[l].append(columns[l][1](data.iloc[j]))
                        j += 1
                    new_data[l] = np.asarray(new_data[l])
                    new_data[len(columns) + num_errs] = np.std(new_data[l], ddof=1) / np.sqrt(np.size(new_data[l]))
                    num_errs = num_errs + 1
                    new_data[l] = (np.asarray(new_data[l])).mean()
                else:
                    new_data[l] = columns[l][1](data.iloc[i])
            destination.loc[k] = new_data
            i = j
            k = k + 1
        self.stats_collections[combined_statistics_collection_name] = destination

    def apply_operation(self, statistics_collection, new_statistics_collection, rules, inputs):
        #[ ] vectorize
        data = self.stats_collections[statistics_collection]
        data_copy = data.copy()
        for i in range(len(data)):
            for j in range(len(rules)):
                data_copy.loc[i, rules[j][0]] = rules[j][1](data_copy.iloc[i], inputs)
        self.stats_collections[new_statistics_collection] = data_copy

    def _calculate_compensation_matrix_n_channels(self, list_comp_samples, channels, threshold=10**-4, k=0.1):
        num_ch = len(channels)
        comp_samples = np.copy(np.asarray(list_comp_samples))
        error_comp_samples = np.copy(comp_samples)
        print(error_comp_samples)

        coeffs = np.diag(np.ones(num_ch))

        errors = np.ones((num_ch, num_ch))

        while not (np.abs(errors) < threshold).all():

            coeffs = np.add(coeffs, np.multiply(-1 * k, errors))

            for i in range(np.shape(coeffs)[0]):
                (error_comp_samples[i])[:, channels] = (np.dot(coeffs, comp_samples[i][:, channels].T)).T

            for i in range(np.shape(coeffs)[0]):
                # error_comp_samples[i][~np.isfinite(error_comp_samples[i])] = 0
                for j in range(np.shape(coeffs)[1]):
                    if i != j:
                        errors[i][j] = (sm.OLS(error_comp_samples[i][:, channels[j]], error_comp_samples[i][:, channels[i]]).fit()).params[0]
                        

        return coeffs

    def _calculate_compensation_matrix_3_channels(self, fcs_ch_1, fcs_ch_2, fcs_ch_3, ch_1, ch_2, ch_3, threshold=10**-4, k=0.1):
        c_12 = c_13 = c_21 = c_23 = c_31 = c_32 = 0.0
        A = np.array([[1.0, c_12, c_13], [c_21, 1.0, c_23], [c_31, c_32, 1.0]])

        copy_fcs_ch_1 = fcs_ch_1.copy()
        copy_fcs_ch_2 = fcs_ch_2.copy()
        copy_fcs_ch_3 = fcs_ch_3.copy()

        e_21 = (sm.OLS(copy_fcs_ch_1[:, ch_2], copy_fcs_ch_1[:, ch_1]).fit()).params[0]
        e_31 = (sm.OLS(copy_fcs_ch_1[:, ch_3], copy_fcs_ch_1[:, ch_1]).fit()).params[0]
        e_12 = (sm.OLS(copy_fcs_ch_2[:, ch_1], copy_fcs_ch_2[:, ch_2]).fit()).params[0]
        e_32 = (sm.OLS(copy_fcs_ch_2[:, ch_3], copy_fcs_ch_2[:, ch_2]).fit()).params[0]
        e_13 = (sm.OLS(copy_fcs_ch_3[:, ch_1], copy_fcs_ch_3[:, ch_3]).fit()).params[0]
        e_23 = (sm.OLS(copy_fcs_ch_3[:, ch_2], copy_fcs_ch_3[:, ch_3]).fit()).params[0]

        n = -1
        if self.full_output:
            fc.plot.density2d(copy_fcs_ch_1, channels=[ch_1, ch_2], mode="scatter")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_1, channels=[ch_1, ch_3], mode="scatter")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_2, channels=[ch_2, ch_1], mode="scatter")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_2, channels=[ch_2, ch_3], mode="scatter")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_3, channels=[ch_3, ch_1], mode="scatter")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_3, channels=[ch_3, ch_2], mode="scatter")
            plt.show()

        while not (np.abs([e_21, e_31, e_12, e_32, e_13, e_23]) < threshold).all():
            if self.full_output:
                n = n + 1
                print("Iteration: " + str(n))
                print("C:")
                print(A)
                print(
                    "\nErrors: "
                    + str(e_12)
                    + ", "
                    + str(e_13)
                    + ", "
                    + str(e_21)
                    + ", "
                    + str(e_23)
                    + ", "
                    + str(e_31)
                    + ", "
                    + str(e_32)
                )
                print("\n\n")
                
            c_21, c_31, c_12, c_32, c_13, c_23 = np.add(
                [c_21, c_31, c_12, c_32, c_13, c_23],
                np.multiply(-1 * k, [e_21, e_31, e_12, e_32, e_13, e_23]),
            )

            A = np.array([[1.0, c_12, c_13], [c_21, 1.0, c_23], [c_31, c_32, 1.0]])

            copy_fcs_ch_1[:, ch_1] = np.dot(
                A[0, :],
                np.asarray(
                    [fcs_ch_1[:, ch_1].T, fcs_ch_1[:, ch_2].T, fcs_ch_1[:, ch_3].T]
                ),
            )
            copy_fcs_ch_1[:, ch_2] = np.dot(
                A[1, :],
                np.asarray(
                    [fcs_ch_1[:, ch_1].T, fcs_ch_1[:, ch_2].T, fcs_ch_1[:, ch_3].T]
                ),
            )
            copy_fcs_ch_1[:, ch_3] = np.dot(
                A[2, :],
                np.asarray(
                    [fcs_ch_1[:, ch_1].T, fcs_ch_1[:, ch_2].T, fcs_ch_1[:, ch_3].T]
                ),
            )
            copy_fcs_ch_2[:, ch_1] = np.dot(
                A[0, :],
                np.asarray(
                    [fcs_ch_2[:, ch_1].T, fcs_ch_2[:, ch_2].T, fcs_ch_2[:, ch_3].T]
                ),
            )
            copy_fcs_ch_2[:, ch_2] = np.dot(
                A[1, :],
                np.asarray(
                    [fcs_ch_2[:, ch_1].T, fcs_ch_2[:, ch_2].T, fcs_ch_2[:, ch_3].T]
                ),
            )
            copy_fcs_ch_2[:, ch_3] = np.dot(
                A[2, :],
                np.asarray(
                    [fcs_ch_2[:, ch_1].T, fcs_ch_2[:, ch_2].T, fcs_ch_2[:, ch_3].T]
                ),
            )
            copy_fcs_ch_3[:, ch_1] = np.dot(
                A[0, :],
                np.asarray(
                    [fcs_ch_3[:, ch_1].T, fcs_ch_3[:, ch_2].T, fcs_ch_3[:, ch_3].T]
                ),
            )
            copy_fcs_ch_3[:, ch_2] = np.dot(
                A[1, :],
                np.asarray(
                    [fcs_ch_3[:, ch_1].T, fcs_ch_3[:, ch_2].T, fcs_ch_3[:, ch_3].T]
                ),
            )
            copy_fcs_ch_3[:, ch_3] = np.dot(
                A[2, :],
                np.asarray(
                    [fcs_ch_3[:, ch_1].T, fcs_ch_3[:, ch_2].T, fcs_ch_3[:, ch_3].T]
                ),
            )

            e_21 = (
                sm.OLS(copy_fcs_ch_1[:, ch_2], copy_fcs_ch_1[:, ch_1]).fit()
            ).params[0]
            e_31 = (
                sm.OLS(copy_fcs_ch_1[:, ch_3], copy_fcs_ch_1[:, ch_1]).fit()
            ).params[0]
            e_12 = (
                sm.OLS(copy_fcs_ch_2[:, ch_1], copy_fcs_ch_2[:, ch_2]).fit()
            ).params[0]
            e_32 = (
                sm.OLS(copy_fcs_ch_2[:, ch_3], copy_fcs_ch_2[:, ch_2]).fit()
            ).params[0]
            e_13 = (
                sm.OLS(copy_fcs_ch_3[:, ch_1], copy_fcs_ch_3[:, ch_3]).fit()
            ).params[0]
            e_23 = (
                sm.OLS(copy_fcs_ch_3[:, ch_2], copy_fcs_ch_3[:, ch_3]).fit()
            ).params[0]

        c_21, c_31, c_12, c_32, c_13, c_23 = np.add(
            [c_21, c_31, c_12, c_32, c_13, c_23],
            np.multiply(-1 * k, [e_21, e_31, e_12, e_32, e_13, e_23]),
        )

        A = np.array([[1.0, c_12, c_13], [c_21, 1.0, c_23], [c_31, c_32, 1.0]])

        if self.full_output:
            n = n + 1
            print("Iteration: " + str(n))
            print("C:")
            print(A)
            print(
                "\nErrors: "
                + str(e_12)
                + ", "
                + str(e_13)
                + ", "
                + str(e_21)
                + ", "
                + str(e_23)
                + ", "
                + str(e_31)
                + ", "
                + str(e_32)
            )
            print("\n\n")
            
            fc.plot.density2d(copy_fcs_ch_1, channels=[ch_1, ch_2], mode="scatter")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_1, channels=[ch_1, ch_3], mode="scatter")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_2, channels=[ch_2, ch_1], mode="scatter")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_2, channels=[ch_2, ch_3], mode="scatter")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_3, channels=[ch_3, ch_1], mode="scatter")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_3, channels=[ch_3, ch_2], mode="scatter")
            plt.show()
            print(A)
        return A

    def _calculate_compensation_matrix_2_channels(self, fcs_ch_1, fcs_ch_2, ch_1, ch_2, threshold=10**-4, k=0.1):
        c_12 = c_21 = 0.0
        A = np.array([[1.0, c_12], [c_21, 1.0]])

        copy_fcs_ch_1 = fcs_ch_1.copy()
        copy_fcs_ch_2 = fcs_ch_2.copy()

        e_21 = (sm.OLS(copy_fcs_ch_1[:, ch_2], copy_fcs_ch_1[:, ch_1]).fit()).params[0]
        e_12 = (sm.OLS(copy_fcs_ch_2[:, ch_1], copy_fcs_ch_2[:, ch_2]).fit()).params[0]

        n = -1
        if self.full_output:
            fc.plot.density2d(copy_fcs_ch_1, channels=[ch_1, ch_2], mode="scatter")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_2, channels=[ch_2, ch_1], mode="scatter")
            plt.show()

        while not (np.abs([e_21, e_12]) < threshold).all():
            if self.full_output:
                n = n + 1
                print("Iteration: " + str(n))
                print("C:")
                print(A)
                print("\nErrors: " + str(e_12) + ", " + str(e_21))
                print("\n\n")
            c_21, c_12 = np.add([c_21, c_12], np.multiply(-1 * k, [e_21, e_12]))

            A = np.array([[1.0, c_12], [c_21, 1.0]])

            copy_fcs_ch_1[:, ch_1] = np.dot(
                A[0, :], np.asarray([fcs_ch_1[:, ch_1].T, fcs_ch_1[:, ch_2].T])
            )
            copy_fcs_ch_1[:, ch_2] = np.dot(
                A[1, :], np.asarray([fcs_ch_1[:, ch_1].T, fcs_ch_1[:, ch_2].T])
            )
            copy_fcs_ch_2[:, ch_1] = np.dot(
                A[0, :], np.asarray([fcs_ch_2[:, ch_1].T, fcs_ch_2[:, ch_2].T])
            )
            copy_fcs_ch_2[:, ch_2] = np.dot(
                A[1, :], np.asarray([fcs_ch_2[:, ch_1].T, fcs_ch_2[:, ch_2].T])
            )

            e_21 = (
                sm.OLS(copy_fcs_ch_1[:, ch_2], copy_fcs_ch_1[:, ch_1]).fit()
            ).params[0]
            e_12 = (
                sm.OLS(copy_fcs_ch_2[:, ch_1], copy_fcs_ch_2[:, ch_2]).fit()
            ).params[0]
        
        (
            c_21,
            c_12,
        ) = np.add([c_21, c_12], np.multiply(-1 * k, [e_21, e_12]))

        A = np.array([[1.0, c_12], [c_21, 1.0]])

        if self.full_output:
            n = n + 1
            print("Iteration: " + str(n))
            print("C:")
            print(A)
            print("\nErrors: " + str(e_12) + ", " + str(e_21))
            print("\n\n")

            fc.plot.density2d(copy_fcs_ch_1, channels=[ch_1, ch_2], mode="scatter")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_2, channels=[ch_2, ch_1], mode="scatter")
            plt.show()
            print(A)
        return A

    def calculate_compensation_matrix(self, sample_collection, compensation_samples, compensation_channels, threshold=10**-4, k=0.1):
        samples_to_compensate = []
        for sample in compensation_samples:
            samples_to_compensate.append(self.sample_collections[sample_collection][sample])
        if len(compensation_channels) == 2:
            self.compensation_matrix = (compensation_channels, self._calculate_compensation_matrix_2_channels(samples_to_compensate[0], samples_to_compensate[1], compensation_channels[0], compensation_channels[1], threshold, k))
        elif len(compensation_channels) == 3:
        #     self.compensation_matrix = self._calculate_compensation_matrix_n_channels(samples_to_compensate[0], samples_to_compensate[1], samples_to_compensate[2], compensation_channels[0], compensation_channels[1], compensation_channels[2], threshold, k)
        # else:
            self.compensation_matrix = (compensation_channels, self._calculate_compensation_matrix_n_channels(samples_to_compensate, compensation_channels, threshold, k))

    def _apply_compensation_matrix_2_channels(self, data_to_compensate, ch_1, ch_2, A):
        data_copy = data_to_compensate.copy()
        for key in data_to_compensate.keys():
            curr_data = data_to_compensate[key]

            curr_data[:, ch_1] = np.dot(
                A[0, :], np.asarray([curr_data[:, ch_1].T, curr_data[:, ch_2].T])
            )
            curr_data[:, ch_2] = np.dot(
                A[1, :], np.asarray([curr_data[:, ch_1].T, curr_data[:, ch_2].T])
            )

            data_copy[key] = curr_data

        return data_copy

    def _apply_compensation_matrix_3_channels(self, data_to_compensate, ch_1, ch_2, ch_3, A):
        data_copy = data_to_compensate.copy()
        for key in data_to_compensate.keys():
            curr_data = data_to_compensate[key]
            curr_data[:, ch_1] = np.dot(A[0, :], np.asarray([curr_data[:, ch_1].T, curr_data[:, ch_2].T, curr_data[:, ch_3].T]))
            curr_data[:, ch_2] = np.dot(A[1, :], np.asarray([curr_data[:, ch_1].T, curr_data[:, ch_2].T, curr_data[:, ch_3].T]))
            curr_data[:, ch_3] = np.dot(A[2, :], np.asarray([curr_data[:, ch_1].T, curr_data[:, ch_2].T, curr_data[:, ch_3].T]))
            data_copy[key] = curr_data
        return data_copy

    def _apply_compensation_matrix_n_channels(self, data_to_compensate, channels, A):
        # [ ] converts channels to tuple
        # [ ] use advanced slicing data[:, (chnls)] = A * data[:, (chnls)]
        data_copy = data_to_compensate.copy()
        for key in data_to_compensate.keys():
            for i in range(len(channels)):
                data_copy[key][i][:, channels] = (np.dot(A, data_copy[:, channels].T)).T
        return data_copy

    def apply_compensation_matrix(self, sample_collection, new_sample_collection):
        compensation_channels = self.compensation_matrix[0]
        compensation_matrix = self.compensation_matrix[1]
        if len(compensation_channels) == 2:
            self.sample_collections[new_sample_collection] = self._apply_compensation_matrix_2_channels(self.sample_collections[sample_collection], compensation_channels[0], compensation_channels[1], compensation_matrix)
        elif len(compensation_channels) == 3:
            self.sample_collections[new_sample_collection] = self._apply_compensation_matrix_3_channels(self.sample_collections[sample_collection], compensation_channels[0], compensation_channels[1], compensation_channels[2], compensation_matrix)
        else:
            self.sample_collections[new_sample_collection] = self._apply_compensation_matrix_n_channels(self.sample_collections[sample_collection], compensation_channels, compensation_matrix)

    def graph_statistics(self, data, errors=(False, False), legend=None, title=None, labels=(None, None), xlog=False, ylog=False, save=True):
        # [ ]: change save to Union[bool, str] so save can be path to file to save
        graphable_data = []
        for _, val in enumerate(data):
            if len(val) <= 3:
                graphable_data.append([self.stats_collections[val[0]], val[1], val[2]])
            else:
                curr_data = self.stats_collections[val[0]]
                for i in range(len(val) - 3):
                    curr_data = curr_data.loc[curr_data[val[3 + i][0]] == val[3 + i][1]]
                graphable_data.append([curr_data, val[1], val[2]])
        plt.close()
        for i in range(len(graphable_data)):
            plt.scatter(graphable_data[i][0][graphable_data[i][1]], graphable_data[i][0][graphable_data[i][2]], zorder=3)
            graph_errors = [None, None]
            for j in range(len(errors)):
                if errors[j] == True:
                    graph_errors[j] = graphable_data[i][0][""+graphable_data[i][1+j]+"_stdErr"]
            plt.errorbar(graphable_data[i][0][graphable_data[i][1]], graphable_data[i][0][graphable_data[i][2]], xerr=graph_errors[0], yerr=graph_errors[1], fmt=" ", capsize = 7, ecolor="#343434", lw=1)
        if xlog:
            plt.xscale('log')
        if ylog:
            plt.yscale('log')
        if legend is not None:
            plt.legend(legend, loc='center left', bbox_to_anchor=(1,0.5))
        plt.title(title, y=1.08)
        plt.xlabel(labels[0])
        plt.ylabel(labels[1])
        if save:
            if title == None:
                title = "untitled"
            title = title.replace('\n', '').replace(':', '-')
            plt.savefig(""+('_').join(('').join(title.split('.')).split(' '))+".png", dpi=500, bbox_inches ="tight")
        plt.show()
   
    def _apply_gate(self, data_to_gate, gating_function, inputs, gate_type):
        data_copy = data_to_gate.copy()
        if 'limits' not in inputs:
            inputs['limits'] = self.lims
        # [ ] all gates to work the same way
        if gate_type != 1:
            return gating_function(data_copy.copy(), r_ready = self.r_ready, **inputs)
        for key in data_copy.keys():
            data_copy[key] = gating_function(data_copy[key], r_ready = self.r_ready, **inputs)
        return data_copy

    def apply_gate(self, sample_collection_to_gate, new_sample_collection, gating_function, inputs = {}, gate_type = 1):
        self.sample_collections[new_sample_collection] = self._apply_gate(self.sample_collections[sample_collection_to_gate], gating_function, inputs, gate_type)
    
    def visualize_plot_change(self, sample_collection_0, data_0, sample_collection_f, data_f, channels: list[str]) -> None:
        self.visualize_plot_overlay([[sample_collection_0, data_0], [sample_collection_f, data_f]], ["#FF0000", "#1E90FF"], channels, [0.04, 0.06], ['.', 'o'])

    def visualize_plot_overlay(self, plots: list[list], colors: list[str], channels: list[str], sizes: list[int] = None, markers: list[str] = None) -> None:
        if sizes is None:
            sizes = [0.01] * len(plots)
        if markers is None:
            markers = ["."] * len(plots)
        for count, plot in enumerate(plots):
            plt.scatter(self.sample_collections[plot[0]][plot[1]][:, channels[0]], self.sample_collections[plot[0]][plot[1]][:, channels[1]], c=colors[count], s=sizes[count], marker=markers[count])
        if channels[0] in self.lims:
            plt.xlim(self.lims[channels[0]])
        if channels[1] in self.lims:
            plt.ylim(self.lims[channels[1]])
        plt.show()
    
