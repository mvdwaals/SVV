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
                        # Trim measurements from here

trim_not_si = np.array([[7540,	160,	2.0,    417,	456,	686,	4.8],
                        [7750,	150,	1.2,    415,	454,	734,	5.5],
                        [8090,	139,	-0.6,   411,	450,	766,	6.6],
                        [8500,	130,	-1.8,   406,	444,	790,	7.7],
                        [7730,	172,	1.5,    424,	463,	817,	4.0],
                        [7290,	180,	3.5,    427,	467,	837,	3.5],
                        [6600,	190,	5.0,    437,	477,	862,	3.0],
                        [6800,	160,	3.0,    430,	470,	897,	4.6]])
