%% Two-Factor learning rule
% (cross-homeo + hybrid homeo/anti-homeo)

function dWdt = kernel_twoTermHybrid(t,W,f_up,params)

	E_set = params.E_set;
	I_set = params.I_set;
	Theta_E = params.Theta_E;
	Theta_I = params.Theta_I;
	alpha = params.alpha;
	beta_E = params.beta_E;
	beta_I = params.beta_I;
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

	dWEEdt = alpha*E*(I_set - I) + beta_E*E*(E_set - E);
	if (W_EE<=0 && dWEEdt<=0)	% if WEE<0 stop updating, unless dWEEdt>0
		dWEEdt = 0;
	end
	dWEIdt = -alpha*I*(I_set - I) - beta_E*I*(E_set - E);
	if (W_EI<=0 && dWEIdt<=0)	% if WEI<0 stop updating, unless dWEIdt>0
		dWEIdt = 0;
	end
	dWIEdt = -alpha*E*(E_set - E) - beta_I*E*(I_set - I);
	if (W_IE<=0 && dWIEdt<=0)	% if WIE<0 stop updating, unless dWIEdt>0
		dWIEdt = 0;
	end
	dWIIdt = alpha*I*(E_set - E) + beta_I*I*(I_set - I);
	if (W_II<=0 && dWIIdt<=0)	% if WII<0 stop updating, unless dWIIdt>0
		dWIIdt = 0;
	end
	dWdt = [dWEEdt,dWEIdt,dWIEdt,dWIIdt];
end


%%
