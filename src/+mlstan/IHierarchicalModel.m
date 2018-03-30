classdef (Abstract) IHierarchicalModel 
	%% IHIERARCHICALMODEL
    %  See also:  M. J. Betancourt and M. Girolami, ArXiv:1312.0906 [Stat] (2013);
    %  S. Weber, A. Gelman, B. Carpenter, D. Lee, M. Betancourt, A. Vehtari, and A. Racine, ArXiv:1602.02055 [Stat] (2016).
    
	%  $Revision$
 	%  was created 08-Dec-2017 16:20:12 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlstan/src/+mlstan.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties (Abstract)
        y        % y_i \sim {\matcal{N}}(\theta_i, \sigma_i^2) measured data
        mathcalD % 1, \ldots, n; data
        theta    % 1, \ldots, n; \theta_i \sim {\mathcal{N}}(\mu, \tau^2); local parameters
        phi      % global parameters
 	end

	methods (Abstract)
        theta = oneWayNormalWithAuxiliary(mu, xi, eta, sigma_eta) 
        % $\theta_i = \mu + \xi\eta_i, \eta_i \sim {\mathcal{N}}(0, \sigma_{\eta}^2)$ with $\tau = |\xi| \sigma_{\eta}$ 
        % and conditioned on auxiliary $\eta$.  However, multiplicative dependence on $\eta$ introduces strong 
        % correlations into joint distribution.      
        
        pi_  = jointDensity_p_q(this) % $\pi(p, q) = \pi(p|q) \pi(q)$
        H_   = Hamiltonian_p_q(this)  % $H(p, q) = -\log\pi(q, q) = -\log\pi(p|q) -\log\pi(q) = T(p|q) + V(q)$
        qdot = dq_dt(this)
        pdot = dp_dt(this)
        
        DV = DeltaV(this) % Hamiltonian transitions
        DT = DeltaT(this) % "
        
        l_ = lambda(this) % to avoid diverging solutions while exploring distributions, step-size must be tuned to curvature
               % $\epsilon \sqrt{\lambda_i} < 2$ for eigenvalues $lambda_i$ of
        M_ = M(this) % $M_{ij} = (\Sigma^{-1})_{ik} \frac{\partial^2 V}{\partial q_k \partial q_j}$
        ch = auxChainAtSmallStepsize(this) % for diagnosing length-scale pathologies
        
        T_ = RiemannT(this) % $T(p,q) = \frac{1}{2}p^T\Sigma^{-1}(q)p - \frac{1}{2}\log|\Sigma(q)|$, metric $\Sigma$
        S_ = Sigma_q(this)  % Riemannian metric ~ Hessian of target distribution
        SA = SoftAbs(this)  % SoftAbs transfomration on Hessian
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

