Active Channels. 

We provide ways to view the activt properties of active channels. One can adapte parameters and save the tuned parameter of active channels for neuron model. 

major parameters to tune: 

1. reversal potential.
2. conductance in the voltage clamp mode; and the time constance in V-clamp mode. 
3. we provide several channelstypes, including voltage gated Na or K channels. 
channames={	'sdNa';...% Somato-dendritic Na+ channel (Na_7)       1
		'aNa';...% Axonic Na+ channel (NaA_2)
		'wCa';...% High thereshold Ca++ channel (Ca_W)    2
                'tCa';...% Low thereshold Ca++ channel (CaT_3)    3
                'aKp';...% A type K+ channel, proximal (K_A_11)   4
                'aKd';...% A type K+ channel, distal (K_A_18)     5
                'drK';...% DR type K+ channel (K_DR_2)            6
                'drKa';...% DR type K+ channel, axonic (K_DRA_4)
                'cK';...% C type K+ channel (K_C_1)               7
                'mK';...% M type K+ channel (K_M_4)               8
                'AHPK';...% AHP type K+ channel (K_AHP_Wtn)       9
                'Ih';...%  h type channel (Ih_3)                  10
                };

4. users can define their own channels. 

Active channels are saved in the general channel kinetics form. 

Active parameters: 
	n: ration of openning 
	alpha: opening rate k_{C to O}
		usually dependens on opening potential V_alpha and constant sigma_alpha
	beta: closing rate  k_{O to C}
		usually dependens on opening potential V_beta and constant sigma_beta
	s.t. in equalliblium state,  n_inf = alpha/(alpha + beta)

	h: inavtivation gating. (in many literature, m is the counterpart of n in K channels)
		parameters governing h is given similar to that for n
	
	
Basic parameters: 
	E_r: reversal potential


by YY: chen@biologie.uni-muenchen.de

