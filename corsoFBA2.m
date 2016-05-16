%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   flux = corsoFBA2(model,onstr,constraint,constrainby,costas)
%   
%   corsoFBA minimizes the overall flux through the entire metabolism
%   according to costs assigned to each reaction. The function first
%   optimizes the objective function then constrains the objective to a given
%   percentage. The model is then simulated again while minimizing the costs
%   given. The function has been tested for maximizing and minimizing a given
%   reaction, but will not work with multiple objectives.
%   
%   This method is described in detail in the following publication:
%     Schultz, A., & Qutub, A. A. (2015). Predicting internal cell fluxes at 
%     sub-optimal growth. BMC systems biology, 9(1), 18.
%   
%   INPUTS:
%     model - metabolic reconstruction to be simulated.
%     onstr - string defining whether the objective function should be
%         maximized or optimize. Options are 'max' and 'min'.
%     constraint - value by which objective will be constrained.
%     constraintby - string. Options 'perc', constraints objective value by
%         percentage from optima, or 'val', constraints flux by absolute 
%         flux value. If constraint is 1, for example, and constraintby is
%         'perc', objective will be constrained to be 1% of its optimal
%         value. If constraintby is 'val', the objective value will be
%         constrained to be 1.
%     costas - cost assigned to each reactions. The options are a numeric
%         value (all reactions will have the same cost), a vector of length
%         model.rxns (reaction model.rxns{i} will have a cost of costas(i)
%         both forward and backwards) or a vector of length twice that of
%         model.rxns (model.rxns{i} will have a cost of costas(i) in the
%         forward direction and a cost of costas(i+length(model.rxns)) in the
%         backwards direction.
%   OUTPUTS:
%     flux - Struct of fields:
%         x - flux distribution with same length as model.rxns
%         y - Dual. Same as in optimizeCbModel
%         f - objective value
%         fm - cost assciate with the flux distribution
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function flux = corsoFBA2(model,onstr,constraint,constrainby,costas)

if strcmp(constrainby,'perc')
    constraint = abs(constraint);
end
%determine objective function
flux1 = optimizeCbModel(model,onstr);
if abs(flux1.f) < 1e-6
%     warning('FBA problem infeasible')
    flux.f = [];
    flux.x = zeros(length(model.rxns),1);
    return
end

%relax results to avoid computational error
if strcmp(constrainby,'perc')
    flux1.f = flux1.x(model.c ~= 0)*(constraint/100);   %Bound
elseif strcmp(constrainby,'val')
    if (flux1.f < constraint) && strcmp(onstr,'max')
        error('Objective Flux not attainable')
    elseif (flux1.f > constraint) && strcmp(onstr,'min')
        error('Objective Flux not attainable')
    else
        flux1.f = constraint;   %Bound
    end
else
    error('Invalid Constraint option');
end
%save original model
model1 = model;

%See if cost is of right length
if length(costas) == 1
    costas = ones(length(model.rxns),1);
end
if ~iscolumn(costas)
    costas = costas';
end
if length(costas)==length(model.rxns)
    costas = [costas; costas];
elseif length(costas) ~= 2*length(model.rxns)
    fprintf('Invalid length of costs\n');
    flux = [];
    return
end

%find internal reactions that are actively reversible
orlen = length(model.rxns);
leng = find(model.lb<0 & model.ub>=0); 

%Tailor model
model.S = [model.S -model.S(:,leng); sparse(zeros(1,orlen+length(leng)))];
model.mets{length(model.mets)+1} = 'pseudomet';
model.b = zeros(length(model.mets),1);
model.c = zeros(orlen+length(leng),1);
model.ub = [model.ub; -model.lb(leng)];
model.lb = zeros(orlen+length(leng),1);
model.S(end,:) = [costas(1:orlen); costas(orlen+leng)];
model.rxns = cat(1,model.rxns,strcat(model.rxns(leng),'added'));

%add reaction for pseudomet consumption
model.rxns = cat(1,model.rxns,'EX_pseudomet');
model.ub = [model.ub; 1e20];
model.lb = [model.lb; 0];
temp = zeros(length(model.mets),1);
temp(end) = -1;
model.S = [model.S temp];
model.c = [model.c; 1];

%change bounds on original optimized reaction
t = find(model1.c);
for k = 1:length(t)
    model = changeRxnBounds(model,model1.rxns(t(k)),flux1.f(k),'b');
    if findRxnIDs(model,[model1.rxns{t(k)} 'added']) ~= 0
        model = changeRxnBounds(model,[model1.rxns{t(k)} 'added'],...
            0,'b');
    end
end

%perform FBA
flux2 = optimizeCbModel(model,'min');

flux.x = flux2.x(1:orlen);
flux.x(leng) = flux.x(leng) - flux2.x((orlen+1):(end-1));
flux.x(abs(flux.x) < 1e-8) = 0;

if isfield(flux1,'y')
    flux.y = flux1.y;
end
if isfield(flux1,'f')
    flux.f = flux1.f;
end
if isfield(flux2,'f')
    flux.fm = flux2.f;
end
end