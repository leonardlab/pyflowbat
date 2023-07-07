import os
import FlowCal as fc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statsmodels.api as sm
import warnings
from typing import Union, Optional, Callable

from . import _std_vals


SampleCollection = dict[str, fc.io.FCSData]
StatisticCollection = pd.DataFrame

comments_to_add = """
        Description

        :param param1: this is a first param
        :param param2: this is a second param
        :returns: this is a description of what is returned
        :raises keyError: raises an exception
        """

class StatisticsExtraction:
    """
    A class containing a rule for extracting statistics from a
    sample collection.

    :param sample_collection_name: the name of the sample collection
        in the workspace from which to extract the statistics
    :type sample_collection_name: str
    :param statistics_collection_name: the name of the statistics
        collection to create/extract statistics into
    :type statistics_collection_name: str
    :param include: the list of keywords in the names of samples
        that must be present to extract statistics from that
        sample, all keywords must be present in a sample's name
        for statistics to be extracted from it
    :type include: list[str]
    :param not_include: the list of keywords in the names of
        samples that must not be present to extract statistics
        from that sample, all keywords must not be present in a
        sample's name for statistics to be extracted from it
    :type not_include: list[str]
    """

    def __init__(
            self,
            sample_collection_name: str,
            statistics_collection_name: str,
            include: list[str],
            not_include: list[str]
            ) -> None:
        """
        Constructor method.
        """
        self.sample_collection_name = sample_collection_name
        self.statistics_collection_name = statistics_collection_name
        self.include = include
        self.not_include = not_include

    def _follows_rule(
            self,
            name: str
        ):
        for word in self.include:
            if word not in name:
                return False
        for word in self.not_include:
            if word in name:
                return False
        return True

class Workspace:
    """
    A class describing PyFlowBAT workspaces. Workspaces
    contain all the methods for operating on batches
    of PyFlowBAT data and contain all samples and
    statistics.
    
    :param stylesheet: the stylesheet for plotting
        PyFlowBAT analyzed data,
        defaults to the PyFlowBAT standard stylesheet
    :type stylesheet: dict
    :param lims_file: the path to the CSV file defining
        upper and lower limits for standard PyFlowBAT
        gating functions
    :type lims_file: str
    :param full_output: whether or not to display all
        output of PyFlowBAT operations;
        it is HIGHLY recommended to leave this at the
        default value:
        True,
        defaults to True
    :type full_output: bool
    """

    ###################################
    # CONSTRUCTOR AND WORKSPACE SETUP #
    ###################################

    def __init__(
            self,
            stylesheet: dict = _std_vals.std_pfb_style,
            lims_file: str = "_std",
            full_output: bool = False
        ) -> None:
        """
        Constructor method.
        """
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
        self.r_ready = False

    def _read_lims_file(self, lims_file) -> dict[str, list[int]]:
        lims_data = pd.read_csv(lims_file)
        lims_dict = {}
        for header in list(lims_data.columns):
            lims_dict[header] = list(lims_data[header])
        return lims_dict

    ###################
    # R FUNCTIONALITY #
    ###################

    def _verify_R_installation(self):
        from subprocess import CalledProcessError
        from subprocess import check_call
        try:
            check_call(['which', 'R'])
        except CalledProcessError:
            raise RuntimeError("No R installation could be found")
        return True

    def init_r(
            self,
            check_R_installation: bool = True
        ) -> None:
        """
        Initialized R functionality for R language gating functions.
        NOTE: this method and R gating functionality both require
        the R programming language to be installed for features to
        work properly.
        
        :param check_R_installation: whether or not to check if the R
            language is installed
            it is HIGHLY recommended to leave this at the
            default value:
            True,
            defaults to True
        :type check_R_installation: bool
        """
        if self.r_ready:
            print("R functionality already initialized")
            return
        if check_R_installation:
            self._verify_R_installation()
        else:
            warnings.warn("Using R gating functions without verifying R is installed is not recommended", RuntimeWarning)
        self.r_ready = True

    ######################
    # BEADS CALCULATIONS #
    ######################

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
                plt.title(f"Beads conversion\n{ch[0]} to  {ch[1]}", y=1.08)
                plt.xlabel(f"{ch[0]} expression")
                plt.ylabel(f"{ch[1]} expression")
                plt.show()
        conv_fctrs = {}
        for count, ch in enumerate(beads_fluorescent_channels):
            conv_fctrs[ch[0]] = models[count].params[0]
            conv_fctrs["" + ch[0] +"_stderr"] = models[count].bse[0]
        return conv_fctrs

    def calculate_beads_factors(
            self,
            beads_file_file_path: str, 
            beads_fluorescent_channels: list[tuple[str, str]],
            beads_num_pops: int,
            beads_conversions_file: str = "_std"
        ) -> None:
        """
        Calculates the conversion factors for this workspace
        from a specified beads file.
        
        :param beads_file_file_path: the path to the beads file
            to calculate conversion factors from
        :type beads_file_file_path: str
        :param beads_fluorescent_channels: a list of tuples containing
            the names of the fluorescent channels and the corresponding
            MEFs;
            follows the pattern of:
            [(\"channel 1\", \"MEF 1\"), (\"channel 2\", \"MEF 2\")]
        :type beads_fluorescent_channels: list[tuple[str, str]]
        :param beads_num_pops: the number of beads populations in the
            beads file to calculate conversion factors from
        :type beads_num_pops: int
        :param beads_conversions_file: the path to the file
            specifying the manufacturers bead conversion values,
            defaults to the PyFlowBAT standard beads values aka
            Spherotech RCP-30-5 Rainbow Calibration Beads
        :type beads_conversions_file: str"""
        self.conversion_factors = self._perform_beads_calculations(
            beads_file_file_path, beads_fluorescent_channels,
            beads_num_pops, beads_conversions_file)

    ##################
    # SAMPLE LOADING #
    ##################

    def load_samples(
            self,
            sample_collection_name: str,
            samples_folder_path: str,
            include: list[str],
            not_include: list[str]
        ) -> None:
        """
        Loads samples from a folder into a sample collection in the workspace.
        
        :param sample_collection_name: the name of the sample collection
            into which to load samples
        :type sample_collection_name: str
        :param samples_folder_path: the path to the folder with the samples
            to load
        :type samples_folder_path: str
        :param include: the list of keywords in the names of files
            that must be present to load, all keywords must be
            present in a file's name to load as a sample
        :type include: list[str]
        :param not_include: the list of keywords in the names of files
            that must not be present to load, all keywords must not be
            present in a file's name to load as a sample
        :type not_include: list[str]
        """
        def samples_quals(name):
            for word in include:
                if word not in name:
                    return False
            for word in not_include:
                if word in name:
                    return False
            return True
        
        extracted_data = {}
        file_names = os.listdir(samples_folder_path)
        file_names.sort()
        for i in range(len(file_names)):
            file_name = str(file_names[i])
            if samples_quals(file_name):
                file_path = os.path.join(samples_folder_path, file_names[i])
                fcs_data = fc.io.FCSData(file_path)
                extracted_data[file_names[i]] = fcs_data 
        self.sample_collections[sample_collection_name] = extracted_data

    #####################################
    # STATISTIC EXTRACTION FROM SAMPLES #
    #####################################
    
    def create_statistic_extraction(
            self,
            sample_collection_name: str,
            statistics_collection_name: str,
            include: list[str],
            not_include: list[str],
            statistic_names: list[str]
    ) -> StatisticsExtraction:
        """
        Creates a rule needed for extracting statistics from a sample collection

        :param sample_collection_name: the name of the sample collection from which to extract
        :type sample_collection_name: str
        :param statistics_collection_name: the name of the statistics collection to create and extract statistics to
        :type statistics_collection_name: str
        :param include: a list of words that must be included in the samples from which to extract statistics
        :type include: list[str]
        :param not_include: a list of words that must NOT be included in the samples from which to extract statistics
        :type not_include: list[str]
        :param statistic_names: a list of the stastics that will be extracted using this rule
        :type statistic_names: list[str]
        :returns: the extraction rule to be used in `extract_samples`
        :rtype: pyflowbat.pyflowbat.StatisticsExtraction
        """
        extraction = StatisticsExtraction(
            sample_collection_name,
            statistics_collection_name,
            include,
            not_include
        )
        df = pd.DataFrame(columns=statistic_names)
        self.stats_collections[statistics_collection_name] = df
        return extraction

    def extract_statistic(
            self,
            extraction: StatisticsExtraction,
            statistc_name: str,
            operation: Callable,
            **kwargs
        ) -> None:
        data = self.sample_collections[extraction.sample_collection_name].copy()
        data_names = list(data.keys())
        data_list = []
        for i in range(len(data_names)):
            file_name = str(data_names[i])
            if extraction._follows_rule(file_name):
                fcs_data = data[file_name]
                data_list.append(operation(name = file_name, data = fcs_data, **kwargs))
        self.stats_collections[extraction.statistics_collection_name][statistc_name] = data_list

    ##########################
    # STATISTIC MANIPULATION #
    ##########################

    def combine_replicates(
            self,
            statistics_collection_name: str,
            combined_statistics_collection_name: str,
            combine_by: list[str],
            combination_operations: Union[str, Callable],
            sem_cols: list[str]
        ) -> None:
        df = self.stats_collections[statistics_collection_name].copy()
        columns_to_drop_combine = [i for i in list(df.columns)
                           if (i not in list(combination_operations.keys())
                               and i not in combine_by)]
        combined = df.drop(
            columns = columns_to_drop_combine
        ).groupby(
            combine_by
        ).aggregate(
            combination_operations
        ).reset_index()
        columns_to_drop_sem = [i for i in list(df.columns)
                           if (i not in sem_cols
                               and i not in combine_by)]
        sems = df.drop(
            columns = columns_to_drop_sem
        ).groupby(
            combine_by
        ).sem(
        ).reset_index(
        ).rename(
            columns = {
                str(x): (str(x) + '_stdErr') for
                x in sem_cols
            }
        )
        df = combined.merge(sems)
        self.stats_collections[combined_statistics_collection_name] = df

    def apply_operation(
            self,
            statistics_collection_name: str,
            new_statistics_collection_name: str,
            statistic_name: str,
            new_statistic_name: str,
            operation: Callable,
            **kwargs
        ) -> None:
        """
        Applies a specified operation to a statistic in a statistics collection.
        
        :param statistics_collection_name: the name of the statistics collection to
            operate on
        :type statistics_collection_name: str
        :param new_statistics_collection_name: the name of the new statistics collection
            to create for the post-operation statistics
        :type new_statistics_collection_name: str
        :param statistic_name: the name of the statistic to operate on
        :type statistic_name: str
        :param new_statistic_name: the name of the statistic to be created from
            the operation
        :type new_statistic_name: str
        :param operation: the function defining the opperation to apply
        :type operation: function/Callable
        :param \*\*kwargs: keywords to pass to the operation
        """
        # note: function can be vectorized
        data = self.stats_collections[statistics_collection_name]
        if new_statistics_collection_name not in self.stats_collections.keys() or self.stats_collections[
            new_statistics_collection_name] is None:
            self.stats_collections[new_statistics_collection_name] = data.copy()
        new_data = self.stats_collections[new_statistics_collection_name]
        if new_statistic_name is None:
            warnings.warn("No new statistic created", RuntimeWarning)
            return
        new_data[new_statistic_name] = (data[statistic_name]).apply(
            operation, axis=1, **kwargs)

    ################
    # COMPENSATION #
    ################

    # def _calculate_compensation_matrix_n_channels(self, list_comp_samples, channels, threshold=10**-4, k=0.1):
    #     num_ch = len(channels)
    #     comp_samples = np.copy(np.asarray(list_comp_samples))
    #     error_comp_samples = np.copy(comp_samples)
    #     print(error_comp_samples)

    #     coeffs = np.diag(np.ones(num_ch))

    #     errors = np.ones((num_ch, num_ch))

    #     while not (np.abs(errors) < threshold).all():

    #         coeffs = np.add(coeffs, np.multiply(-1 * k, errors))

    #         for i in range(np.shape(coeffs)[0]):
    #             (error_comp_samples[i])[:, channels] = (np.dot(coeffs, comp_samples[i][:, channels].T)).T

    #         for i in range(np.shape(coeffs)[0]):
    #             # error_comp_samples[i][~np.isfinite(error_comp_samples[i])] = 0
    #             for j in range(np.shape(coeffs)[1]):
    #                 if i != j:
    #                     errors[i][j] = (sm.OLS(error_comp_samples[i][:, channels[j]], error_comp_samples[i][:, channels[i]]).fit()).params[0]
                        

    #     return coeffs

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
            plt.title(f"Pre-compensation scatter: impact of\n{ch_1} on {ch_2}")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_1, channels=[ch_1, ch_3], mode="scatter")
            plt.title(f"Pre-compensation scatter: impact of\n{ch_1} on {ch_3}")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_2, channels=[ch_2, ch_1], mode="scatter")
            plt.title(f"Pre-compensation scatter: impact of\n{ch_2} on {ch_2}")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_2, channels=[ch_2, ch_3], mode="scatter")
            plt.title(f"Pre-compensation scatter: impact of\n{ch_2} on {ch_3}")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_3, channels=[ch_3, ch_1], mode="scatter")
            plt.title(f"Pre-compensation scatter: impact of\n{ch_3} on {ch_1}")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_3, channels=[ch_3, ch_2], mode="scatter")
            plt.title(f"Pre-compensation scatter: impact of\n{ch_3} on {ch_2}")
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
            plt.title(f"Post-compensation scatter: impact of\n{ch_1} on {ch_2}")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_1, channels=[ch_1, ch_3], mode="scatter")
            plt.title(f"Post-compensation scatter: impact of\n{ch_1} on {ch_3}")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_2, channels=[ch_2, ch_1], mode="scatter")
            plt.title(f"Post-compensation scatter: impact of\n{ch_2} on {ch_1}")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_2, channels=[ch_2, ch_3], mode="scatter")
            plt.title(f"Post-compensation scatter: impact of\n{ch_2} on {ch_3}")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_3, channels=[ch_3, ch_1], mode="scatter")
            plt.title(f"Post-compensation scatter: impact of\n{ch_3} on {ch_1}")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_3, channels=[ch_3, ch_2], mode="scatter")
            plt.title(f"Post-compensation scatter: impact of\n{ch_3} on {ch_2}")
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
            plt.title(f"Pre-compensation scatter: impact of\n{ch_1} on {ch_2}")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_2, channels=[ch_2, ch_1], mode="scatter")
            plt.title(f"Pre-compensation scatter: impact of\n{ch_2} on {ch_1}")
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
            plt.title(f"Post-compensation scatter: impact of\n{ch_1} on {ch_2}")
            plt.show()
            fc.plot.density2d(copy_fcs_ch_2, channels=[ch_2, ch_1], mode="scatter")
            plt.title(f"Post-compensation scatter: impact of\n{ch_2} on {ch_1}")
            plt.show()
            print(A)
        return A

    def calculate_compensation_matrix(
            self,
            sample_collection_name: str,
            compensation_samples_name: str,
            compensation_channel_names: list[str],
            threshold: int = 10**-4,
            compensation_rate: float = 0.1
        ) -> None:
        """
        Calculates a compensation matrix by repeatedly flattening to zero
        the best fit line through the specified compensation sample data.
        Flattening stops once the best fit line has a slope under the provided
        threshold.
        NOTE: at present, only 2- and 3- sample compensation matrices can
        be calculated.
        NOTE: at present, flattening continues for all compensation samples
        as long as any compensation sample has a best fit slope greater than
        the threshold.
        NOTE: this method must always be run before the workspace's
        `apply_compensation_matrix` method.
        
        :param sample_collection_name: the name of the collection from
            which to find the samples for compensation matrix calculation
        :type sample_collection_name: str
        :param compensation_sample_names: the names of the samples used
            for compensation
        :type compensation_sample_names: list[str]
        :param compensation_channel_names: the channel names to be
            compensated
        :type compensation_channel_names: list[str]
        :param threshold: the threshold to flatten the best fit
            lines to,
            defaults to 10**4
        :type threshold: float
        :param compensation_rate: the rate at which to calculate the
            compensation matrix,
            defaults to 0.1
        :type compensation_rate: float
        """
        samples_to_compensate = []
        for sample in compensation_samples_name:
            samples_to_compensate.append(self.sample_collections[sample_collection_name][sample])
        if len(compensation_channel_names) == 2:
            self.compensation_matrix = (compensation_channel_names, self._calculate_compensation_matrix_2_channels(
                samples_to_compensate[0], samples_to_compensate[1], compensation_channel_names[0], compensation_channel_names[1], threshold, compensation_rate))
        elif len(compensation_channel_names) == 3:
            self.compensation_matrix = self._calculate_compensation_matrix_n_channels(
                samples_to_compensate[0], samples_to_compensate[1], samples_to_compensate[2],
                compensation_channel_names[0], compensation_channel_names[1], compensation_channel_names[2],
                threshold, compensation_rate)
        else:
            raise NotImplementedError("Compensation not implemented for more than 3 colors") 
        #     self.compensation_matrix = (compensation_channels, self._calculate_compensation_matrix_n_channels(samples_to_compensate, compensation_channels, threshold, k))

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

    # def _apply_compensation_matrix_n_channels(self, data_to_compensate, channels, A):
    #     # note: function can be vectorized
    #     data_copy = data_to_compensate.copy()
    #     for key in data_to_compensate.keys():
    #         for i in range(len(channels)):
    #             data_copy[key][i][:, channels] = (np.dot(A, data_copy[:, channels].T)).T
    #     return data_copy

    def apply_compensation_matrix(
            self,
            sample_collection_name: str,
            new_sample_collection_name: str
        ) -> None:
        """
        Applies the workspace compensation matrix to one sample collection creating another.
        NOTE: the workspace's `calculate_compensation_matrix` method must be run before this
        method may be run.
        
        :param sample_collection_name: the name of the sample collection to compensate
        :type sample_collection_name: str
        :param new_sample_collection_name: the name of the new sample collection to create
            for the compensated samples
        :type new_sample_collection_name: str
        """
        if self.compensation_matrix is None:
            raise ValueError("The workspace compensation matrix is not defined, please run calculate a matrix first")
        compensation_channels = self.compensation_matrix[0]
        compensation_matrix = self.compensation_matrix[1]
        if len(compensation_channels) == 2:
            self.sample_collections[new_sample_collection_name] = self._apply_compensation_matrix_2_channels(
                self.sample_collections[sample_collection_name], compensation_channels[0],
                compensation_channels[1], compensation_matrix)
        elif len(compensation_channels) == 3:
            self.sample_collections[new_sample_collection_name] = self._apply_compensation_matrix_3_channels(
                self.sample_collections[sample_collection_name], compensation_channels[0],
                compensation_channels[1], compensation_channels[2], compensation_matrix)
        else:
            raise NotImplementedError("Compensation not implemented for more than 3 colors") 
            # self.sample_collections[new_sample_collection] = self._apply_compensation_matrix_n_channels(self.sample_collections[sample_collection], compensation_channels, compensation_matrix)

    ##########
    # GATING #
    ##########

    def apply_gate(
            self,
            sample_collection_name: str,
            new_sample_collection_name: str,
            gating_function: Callable,
            **kwargs
        ) -> None:
        """
        Applies a specified gate to a sample collection.
        
        :param sample_collection_name: the name of the sample collection to gate
        :type sample_collection_name: str
        :param new_sample_collection_name: the name of the new sample collection to create
            for the gated samples
        :type new_sample_collection_name: str
        :param gating_function: the function defining the gate to apply
        :type gating_function: function/Callable
        :param \*\*kwargs: keywords to pass to the gating function
        """
        data_copy = (self.sample_collections[sample_collection_name]).copy()
        self.sample_collections[new_sample_collection_name] = gating_function(
            data_copy.copy(), r_ready = self.r_ready, limits = self.lims, **kwargs)

    #################
    # VISUALIZATION #
    #################

    def graph_statistics(
            self,
            data: list[list[Union[str, list[str]]]],
            errors: tuple[bool, bool] = (False, False),
            legend: Optional[str] = None,
            title: Optional[str] = None,
            labels: tuple[Optional[str], Optional[str]] = (None, None),
            xlog: bool = False,
            ylog: bool = False,
            save: Union[bool, str] = True
        ) -> None:
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
        if save == True or isinstance(save, str):
            if title == None:
                save_path = "untitled"
            save_path = title.replace('\n', '').replace(':', '-')
            if isinstance(save, str):
                save_path = save
            plt.savefig(""+('_').join(('').join(save_path.split('.')).split(' '))+".png", dpi=500, bbox_inches ="tight")
        plt.show()
    
    def visualize_plot_change(
            self,
            sample_collection_name_0: str,
            sample_name_0: str,
            sample_collection_name_f: str,
            sample_name_f: str,
            channels: list[str]
        ) -> None:
        self.visualize_plot_overlay(
            [[sample_collection_name_0, sample_name_0], [sample_collection_name_f, sample_name_f]],
            ["#FF0000", "#1E90FF"],
            channels,
            [0.04, 0.06],
            ['.', 'o'],
            [sample_collection_name_0, sample_collection_name_f],
            f"Change in sample {sample_name_0}\naka {sample_name_f}: {sample_collection_name_0} to {sample_collection_name_f}")

    def visualize_plot_overlay(
            self,
            plots: list[list],
            colors: list[str],
            channels: list[str],
            sizes: Optional[list[int]] = None,
            markers: Optional[list[str]] = None,
            legend: Optional[list] = None,
            title: Optional[str] = None
        ) -> None:
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
        if legend is not None:
            plt.legend(legend, loc='center left', bbox_to_anchor=(1,0.5))
        if title is not None:
            plt.title(title, y=1.08)
        plt.xlabel(channels[0])
        plt.ylabel(channels[1])
        plt.show()
    
