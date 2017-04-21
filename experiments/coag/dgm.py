from matplotlib import rc
rc("font", family = "serif", size = 10)
rc("text", usetex=True)

import daft

h_shift = 1.5
v_shift = 2

pgm = daft.PGM([15.0, 10.0], origin=[0, 0])

# CCA Latent
pgm.add_node(daft.Node("z_bg", r"$z_{bc}$", 1 + h_shift, 4 + v_shift))
pgm.add_node(daft.Node("z_vitals", r"$z_{vt}$", 3 + h_shift, 4 + v_shift))
pgm.add_node(daft.Node("z_common", r"$z_{cmn}$", 5 + h_shift, 4 + v_shift))
pgm.add_node(daft.Node("z_cbc", r"$z_{cbc}$", 7 + h_shift, 4 + v_shift))
pgm.add_node(daft.Node("z_cascade", r"$z_{csc}$", 9 + h_shift, 4 + v_shift))
pgm.add_node(daft.Node("z_pt", r"$z_{pt}$", 11 + h_shift, 4 + v_shift))

# CCA Shared Latent Matrices
pgm.add_node(daft.Node("W_bg", r"$W_{bg}$", 2 + h_shift, 5.5 + v_shift))
pgm.add_node(daft.Node("W_vitals", r"$W_{vt}$", 4 + h_shift, 5.5 + v_shift))
pgm.add_node(daft.Node("W_cbc", r"$W_{cbc}$", 6 + h_shift, 5.5 + v_shift))
pgm.add_node(daft.Node("W_cascade", r"$W_{csc}$", 8 + h_shift, 5.5 + v_shift))
pgm.add_node(daft.Node("W_pt", r"$W_{pt}$", 10 + h_shift, 5.5 + v_shift))

# CCA Hierarchical Shared Latent Matrices
pgm.add_node(daft.Node("W_0_bg", r"$W_{0,bg}$", 2 + h_shift, 6.5 + v_shift))
pgm.add_node(daft.Node("W_0_vitals", r"$W_{0,vt}$", 4 + h_shift, 6.5 + v_shift))
pgm.add_node(daft.Node("W_0_cbc", r"$W_{0,cbc}$", 6 + h_shift, 6.5 + v_shift))
pgm.add_node(daft.Node("W_0_cascade", r"$W_{0,csc}$", 8 + h_shift, 6.5 + v_shift))
pgm.add_node(daft.Node("W_0_pt", r"$W_{0,pt}$", 10 + h_shift, 6.5 + v_shift))

# CCA Individual Latent Matrices
pgm.add_node(daft.Node("B_bg", r"$B_{bg}$", 2 + h_shift, 0.5 + v_shift))
pgm.add_node(daft.Node("B_vitals", r"$B_{vts}$", 4 + h_shift, 0.5 + v_shift))
pgm.add_node(daft.Node("B_cbc", r"$B_{cbc}$", 6 + h_shift, 0.5 + v_shift))
pgm.add_node(daft.Node("B_cascade", r"$B_{csc}$", 8 + h_shift, 0.5 + v_shift))
pgm.add_node(daft.Node("B_pt", r"$B_{pt}$", 10 + h_shift, 0.5 + v_shift))

pgm.add_node(daft.Node("B_0_bg", r"$B_{0,bg}$", 2 + h_shift, -0.5 + v_shift))
pgm.add_node(daft.Node("B_0_vitals", r"$B_{0,vts}$", 4 + h_shift, -0.5 + v_shift))
pgm.add_node(daft.Node("B_0_cbc", r"$B_{0,cbc}$", 6 + h_shift, -0.5 + v_shift))
pgm.add_node(daft.Node("B_0_cascade", r"$B_{0,csc}$", 8 + h_shift, -0.5 + v_shift))
pgm.add_node(daft.Node("B_0_pt", r"$B_{0,pt}$", 10 + h_shift, -0.5 + v_shift))

# Observed
pgm.add_node(daft.Node("x_bg", r"$x_{bg}$", 2 + h_shift, 2 + v_shift, observed = True))
pgm.add_node(daft.Node("x_vitals", r"$x_{vts}$", 4 + h_shift, 2 + v_shift, observed = True))
pgm.add_node(daft.Node("x_cbc", r"$x_{cbc}$", 6 + h_shift, 2 + v_shift, observed = True))
pgm.add_node(daft.Node("x_cascade", r"$x_{csc}$", 8 + h_shift, 2 + v_shift, observed = True))
pgm.add_node(daft.Node("x_pt", r"$x_{pt}$", 10 + h_shift, 2 + v_shift, observed = True))

# Edges

# Latent to observed
pgm.add_edge("z_common", "x_bg")
pgm.add_edge("z_common", "x_vitals")
pgm.add_edge("z_common", "x_cbc")
pgm.add_edge("z_common", "x_cascade")
pgm.add_edge("z_common", "x_pt")

pgm.add_edge("z_bg", "x_bg")
pgm.add_edge("z_vitals", "x_vitals")
pgm.add_edge("z_cbc", "x_cbc")
pgm.add_edge("z_cascade", "x_cascade")
pgm.add_edge("z_pt", "x_pt")

# Latent Shared to observed
pgm.add_edge("W_bg", "x_bg")
pgm.add_edge("W_vitals", "x_vitals")
pgm.add_edge("W_cbc", "x_cbc")
pgm.add_edge("W_cascade", "x_cascade")
pgm.add_edge("W_pt", "x_pt")

# Prior to Latent Shared
pgm.add_edge("W_0_bg", "W_bg")
pgm.add_edge("W_0_vitals", "W_vitals")
pgm.add_edge("W_0_cbc", "W_cbc")
pgm.add_edge("W_0_cascade", "W_cascade")
pgm.add_edge("W_0_pt", "W_pt")

# Latent Individual to observed
pgm.add_edge("B_bg", "x_bg")
pgm.add_edge("B_vitals", "x_vitals")
pgm.add_edge("B_cbc", "x_cbc")
pgm.add_edge("B_cascade", "x_cascade")
pgm.add_edge("B_pt", "x_pt")

# Prior to Individual Latent
pgm.add_edge("B_0_bg", "B_bg")
pgm.add_edge("B_0_vitals", "B_vitals")
pgm.add_edge("B_0_cbc", "B_cbc")
pgm.add_edge("B_0_cascade", "B_cascade")
pgm.add_edge("B_0_pt", "B_pt")

# Plates
pgm.add_plate(daft.Plate([0 + h_shift, 1 + v_shift, 12, 4], label=r"$N_j$", rect_params = {'joinstyle':'round'}))

# Hierarchical Plate
pgm.add_plate(daft.Plate([-1 + h_shift, 0 + v_shift, 14, 6], label=r"$J$", rect_params = {'joinstyle':'round'}))

pgm.render()
pgm.figure.savefig("CCA.png", dpi = 1000)
