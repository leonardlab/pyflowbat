import rpy2.robjects as ro

r = ro.r
r['source']('rp2_test.r')

r_func = ro.globalenv['n_to_power']

r_result = r_func(4, 5)

print(r_result)