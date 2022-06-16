%% StandardHomeo learning rule


function dWdt = kernel_standardHomeo(t,W,f_up,params)

	E_set = params.E_set;
	I_set = params.I_set;
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
	
	if (W_EE<0)
		W_EE = 0;
	end
	if (W_EI<0)
		W_EI = 0;
	end
	if (W_IE<0)
		W_IE = 0;
	end
	if (W_II<0)
		W_II = 0;
	end
	E_aux = f_up{1}(W_EE,W_EI,W_IE,W_II);
	I_aux = f_up{2}(W_EE,W_EI,W_IE,W_II);
	if (W_EE*E_aux-W_EI*I_aux >= Theta_E)
		E = E_aux;	% suprathreshold
	else
		E = 0;	% subthreshold
	end
	if (W_IE*E_aux-W_II*I_aux >= Theta_I)
		I = I_aux;
	else
		I = 0;
	end
	
	dWEEdt = alpha_EE*E*(E_set - E);
	dWEIdt = -alpha_EI*I*(E_set - E);
	dWIEdt = alpha_IE*E*(I_set - I);
	dWIIdt = -alpha_II*I*(I_set - I);
	dWdt = [dWEEdt,dWEIdt,dWIEdt,dWIIdt];
end


%%
