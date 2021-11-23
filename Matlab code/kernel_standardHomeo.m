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
	E_aux = f_up{1}(W_EE,W_EI,W_IE,W_II);
	I_aux = f_up{2}(W_EE,W_EI,W_IE,W_II);
	if (W_EE*E_aux-W_EI*I_aux >= Theta_E)
		E = E_aux;
	else
		E = 0;
	end
	if (W_IE*E_aux-W_II*I_aux >= Theta_I)
		I = I_aux;
	else
		I = 0;
	end
	
	dWEEdt = alpha_EE*E*(E_set - E);
	if (W_EE<=0 && dWEEdt<=0)	% if WEE<0 stop updating, unless dWEEdt>0
		dWEEdt = 0;
	end
	dWEIdt = -alpha_EI*I*(E_set - E);
	if (W_EI<=0 && dWEIdt<=0)	% if WEI<0 stop updating, unless dWEIdt>0
		dWEIdt = 0;
	end
	dWIEdt = alpha_IE*E*(I_set - I);
	if (W_IE<=0 && dWIEdt<=0)	% if WIE<0 stop updating, unless dWIEdt>0
		dWIEdt = 0;
	end
	dWIIdt = -alpha_II*I*(I_set - I);
	if (W_II<=0 && dWIIdt<=0)	% if WII<0 stop updating, unless dWIIdt>0
		dWIIdt = 0;
	end
	dWdt = [dWEEdt,dWEIdt,dWIEdt,dWIIdt];
end


%%
