%% MADEXbrussodeJacComp: Brusselator Compressed Jacobian
% This example shows how to use the |fmad| class to calculate the Jacobian
% of the Brusselator problem using compressed storage. If a good coloring
% can be found then this approach is usually more efficient than use of
% sparse or full storage as in |MADEXbrussodeJac| and |MADEXbrussodeJacSparse|.
%
% See also: brussode_f, brussode_S, MADEXbrussodeJac, MADEXbrussodeJacSparse, brussode

%% Initialise variables
% We define the input variables, setting N to 10 for this example.
format compact
N=10;
t=0;
y0=ones(2*N,1);

%% Obtain the sparsity pattern
% The function brussode_S returns the sparsity pattern, the matrix with an
% entry 1 when the corresponding Jacobian entry is nonzero and entry 0 when
% the corresponding Jacobian entry is always zero.
sparsity_pattern=brussode_S(N);

%% Determine the coloring
% The function MADcolor takes the sparsity pattern and determines a color
% (or group number) for each component of the input variables so that
% perturbing one component with that color affects different dependent
% variables to those affected by any other variable in that group.
color_groups=dsjmcolor(sparsity_pattern);
%sparsity_pattern
%color_groups
ncolors=max(color_groups) % number of colors used
