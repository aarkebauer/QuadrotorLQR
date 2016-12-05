% Determines whether an input system described as (A,B,C) is controllable and/or observable
% A. Arkebauer, D. Cody
% EENG 436 Final Project
% 12/4/16

function contObs(A,B,C)

%% Controllability
[~, dim] = size(A);

for i=1:dim
    if i == 1
        contMat = B;
    else
        contMat = [contMat A^(i-1)*B];
    end
end

if rank(contMat) == dim % if contMat is full rank
    disp('System is Controllable')
else
    disp(['System is not Controllable. Rank = ' num2str(rank(contMat)) ' (must be ' num2str(dim) ')'])
end


%% Observability
for i=1:dim
    if i == 1
        obsMat = C;
    else
        obsMat = [obsMat; C*A^(i-1)];
    end
end

if rank(obsMat) == dim % if contMat is full rank
    disp('System is Observable')
else
    disp(['System is not Observable. Rank = ' num2str(rank(obsMat)) ' (must be ' num2str(dim) ')'])
end


end