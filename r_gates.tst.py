
import pyflowbat as pfb

wrksp = pfb.Workspace()

wrksp.init_r()

ungated_data = pfb.fc.io.FCSData("./test-files/PreDox_B_002.fcs")
gated_data = pfb.r_gating.singlet_gate(ungated_data, ["FSC-A", "FSC-H"], r_ready = True)
print(gated_data.shape)
gated_data = pfb.r_gating.clust_2d_gate(ungated_data, ["FSC-A", "SSC-A"], target = [7.5*10**4, 5*10**4], r_ready = True)
print(gated_data.shape)
gated_data = pfb.r_gating.transitional_gate(ungated_data, ["FSC-A", "SSC-A"], target = [7.5*10**4, 5*10**4], r_ready = True)
print(gated_data.shape)