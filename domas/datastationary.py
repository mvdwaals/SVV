import numpy as np
#   0       1       2       3       4           5           6
# hp[ft] Vc[kts] Tmta[c] FFl[lbhr] FFr[lbhr] mfu[lb] alpha[deg]
data_not_si = np.array([#[5010, 248, 11.2, 768, 784, 215, 1.5], # This is a test measurement
                        [6040, 251, 10.5, 786, 776, 410, 1.3],
                        [6035, 222, 8.2, 637, 650, 450, 2.0],
                        [6030, 192, 6.1, 515, 548, 487, 3.0], # Thrust here seems off
                        [6040, 165, 4.8, 558, 470, 540, 4.5],
                        [6030, 130, 3.5, 373, 385, 574, 8.3],
                        [6040, 115, 3.0, 414, 435, 604, 10.1]])
