%% CrossHomeo learning rule


function dWdt = kernel_crossHomeo(t,W,f_up,stable_conds,params)

	E_set = params.E_set;
	I_set = params.I_set;
	E_max = params.E_max;
	I_max = params.I_max;
	E_min = params.E_min;
	I_min = params.I_min;
	Theta_E = params.Theta_E;
	Theta_I = params.Theta_I;
	alpha_EE = params.alpha_EE;
	alpha_EI = params.alpha_EI;
	alpha_IE = params.alpha_IE;
	alpha_II = params.alpha_II;
	W_EE = W(1);
	W_EI = W(2);
	W_IE = W(3);
	W_II = W(4);
	
% 	if (W_EE<0)
% 		W_EE = 0;
% 	end
% 	if (W_EI<0)
% 		W_EI = 0;
% 	end
% 	if (W_IE<0)
% 		W_IE = 0;
% 	end
% 	if (W_II<0)
% 		W_II = 0;
% 	end
	E_aux = f_up{1}(W_EE,W_EI,W_IE,W_II);
	I_aux = f_up{2}(W_EE,W_EI,W_IE,W_II);
	Det = stable_conds{1}(W_EE,W_EI,W_IE,W_II);
	Tr = stable_conds{2}(W_EE,W_EI,W_IE,W_II);
% 	E = E_aux;
% 	I = I_aux;
% 	if (W_EE*E_aux-W_EI*I_aux >= Theta_E && W_IE*E_aux-W_II*I_aux >= Theta_I)
	if (W_EE*E_aux-W_EI*I_aux >= Theta_E && W_IE*E_aux-W_II*I_aux >= Theta_I)
		supra_flag = 1;
		if (Det>0 && Tr<0)	% stable
			E = E_aux;	% suprathreshold
			I = I_aux;
		else
			E = E_min;
			I = I_min;
		end
	else
		supra_flag = 0;	% subthreshold
		if (Det<=0)
			E = E_max;
			I = I_max;
		else
			E = E_min;
			I = I_min;
		end
	end
% % 	if (W_IE*E_aux-W_II*I_aux >= Theta_I)
% % 		I = I_aux;
% % 	else
% % 		I = I_min;
% % 	end
% 
	if (Det<=0 || Tr>=0)	% unstable
		if supra_flag==1
			E = E_min;
			I = I_min;
		else
			E = E_max;
			I = I_max;
		end
	end

	dWEEdt = alpha_EE*E*(I_set - I);
	dWEIdt = -alpha_EI*I*(I_set - I);
	dWIEdt = -alpha_IE*E*(E_set - E);
	dWIIdt = alpha_II*I*(E_set - E);
	dWdt = [dWEEdt,dWEIdt,dWIEdt,dWIIdt];

end


%%
