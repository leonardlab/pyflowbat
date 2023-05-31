# import numpy as np

# def convert_to_mefs(
#         # mef_col_name: str,
#         # conversion_factors_name: str
# ):
#     lambda row, inputs: row['Mean FITC-A'] * my_wrkspc.conversion_factors["FITC-A"] # [ ] figure out how to deal with this

# def _convert_to_mefs(
#         # arbitrary_unit_col_name: str,
#         # mef_col_name: str,
#         # conversion_factors_name: str
# ):
#     pass

# def calculate_std_error():
#     pass


# class{
#         vars

#         function(){
#             func2()
#         }
# }

# func2(){
#     vars
# }

# my_wrkspc.apply_operation('dose response combined', 'dose response converted', 
#     [
#         ['', convert_to_mefs, {}]
#         ['MEFLs', lambda row, inputs: row['Mean FITC-A']*my_wrkspc.conversion_factors["FITC-A"]],
#         ['MEFLs_stdErr', lambda row, inputs: row['Mean FITC-A']*my_wrkspc.conversion_factors["FITC-A"] * np.sqrt((row['Mean FITC-A_stdErr']/row['Mean FITC-A'])**2+(my_wrkspc.conversion_factors["FITC-A_stderr"]/my_wrkspc.conversion_factors["FITC-A"])**2)],
#         ['MEBFPs', lambda row, inputs: row['Mean Pacific Blue-A']*my_wrkspc.conversion_factors["Pacific Blue-A"]],
#         ['MEBFPs_stdErr', lambda row, inputs: row['Mean Pacific Blue-A']*my_wrkspc.conversion_factors["Pacific Blue-A"] * np.sqrt((row['Mean Pacific Blue-A_stdErr']/row['Mean Pacific Blue-A'])**2+(my_wrkspc.conversion_factors["Pacific Blue-A_stderr"]/my_wrkspc.conversion_factors["Pacific Blue-A"])**2)],
#         ['MEFLs (sub)', lambda row, inputs: my_wrkspc.conversion_factors["FITC-A"] * (row['Mean FITC-A'] - my_wrkspc.stats_collections["tx control combined"]["Mean FITC-A"][0])],
#         ['MEFLs (sub)_stdErr', lambda row, inputs: my_wrkspc.conversion_factors["FITC-A"] * (row['Mean FITC-A'] - my_wrkspc.stats_collections["tx control combined"]["Mean FITC-A"][0]) * np.sqrt( (my_wrkspc.conversion_factors["FITC-A_stderr"]/my_wrkspc.conversion_factors["FITC-A"])**2 + ( np.sqrt( (row['Mean FITC-A_stdErr'])**2 + (my_wrkspc.stats_collections["tx control combined"]["Mean FITC-A"][0])**2 ) / (row['Mean FITC-A'] - my_wrkspc.stats_collections["tx control combined"]["Mean FITC-A"][0]) )**2 )]
#     ], []
# )