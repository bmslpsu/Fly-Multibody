function [equations] = controller_from_cl(saveB)
%% Inputs
if nargin < 1
    saveB = false;
end

%% Derive the head & body controller expressions from corresponding plants and closed-loop models
syms C1 C2 P1 P2 H1 H2

G1 = C1*P1;
G2 = C2*P2;

eq1 = H1 == G1 / (1 + G1 + G2);
eq2 = H2 == G2 / (1 + G1 + G2);

equations.multi.closedloop.H1 = G1 / (1 + G1 + G2);
equations.multi.closedloop.H2 = G2 / (1 + G1 + G2);
equations.multi.controllers = solve([eq1 eq2], [C1 C2]);

%% Derive the body controller expression from corresponding plant and closed-loop model
syms C1 C2 P1 P2 H1 H2

G1 = C1*P1;

eq1 = H1 == G1 / (1 + G1);

equations.single.closedloop.H1 = G1 / (1 + G1);
C1_single = solve(eq1, C1);
equations.single.controllers.C1 = C1_single;

%% Save
if saveB
    fname = 'Symbolic_expressions'; 
    root = 'E:\DATA\Magno_Data\Multibody';
    savedir = fullfile(root,'sym');
    mkdir(savedir)
    save(fullfile(savedir, [fname '.mat']), 'equations');
end

end