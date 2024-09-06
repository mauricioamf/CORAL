function rxnsToAdd = copyRxnInfo(ecModel, rxnsToAdd, index)
% copyRxnInfo
%
%   Copy annotations from an old reaction into a new reaction. This 
%   function is not intended to be called outside of 'expandEnzComplex' or 
%   for any other usage. Code adapted from the 'changeRxns' function from
%   the RAVEN Toolbox.
%
%   Usage
%       rxnsToAdd = copyRxnInfo(ecModel, rxnsToAdd, index)
%
%   Parameters
%       ecModel                       a GECKO3 ecModel structure
%       rxnsToAdd                     struct containing rxns, mets and
%                                     stoichCoeffs for reactions being
%                                     altered
%       index                         index of original reaction
%
%   Outputs
%       rxnsToAdd                     struct containing rxns, mets and
%                                     stoichCoeffs for reactions being
%                                     altered, now with annotations such as
%                                     eccodes, miriams, and others.
%
% .. Author:
%       - Mauricio Ferreira &         2023.09.29
%         Eduardo Almeida

if isfield(ecModel,'rxnNames')
    rxnsToAdd.rxnNames=ecModel.rxnNames(index);
end
if isfield(ecModel,'lb')
    rxnsToAdd.lb=ecModel.lb(index);
end
if isfield(ecModel,'ub')
    rxnsToAdd.ub=ecModel.ub(index);
end
if isfield(ecModel,'c')
    rxnsToAdd.c=ecModel.c(index);
end
if isfield(ecModel,'eccodes')
    rxnsToAdd.eccodes=ecModel.eccodes(index);
end
if isfield(ecModel,'subSystems')
    rxnsToAdd.subSystems=ecModel.subSystems(index);
end
if isfield(ecModel,'rxnComps')
    rxnsToAdd.rxnComps=ecModel.rxnComps(index);
end
% if isfield(model,'grRules')
%     rxnsToAdd.grRules=model.grRules(index);
% end
if isfield(ecModel,'rxnFrom')
    rxnsToAdd.rxnFrom=ecModel.rxnFrom(index);
end
if isfield(ecModel,'rxnScores')
    rxnsToAdd.rxnScores=ecModel.rxnScores(index);
end
if isfield(ecModel,'rxnMiriams')
    rxnsToAdd.rxnMiriams=ecModel.rxnMiriams(index);
end
if isfield(ecModel,'rxnNotes')
    rxnsToAdd.rxnNotes=ecModel.rxnNotes(index);
end
if isfield(ecModel,'rxnReferences')
    rxnsToAdd.rxnReferences=ecModel.rxnReferences(index);
end
if isfield(ecModel,'rxnConfidenceScores')
    rxnsToAdd.rxnConfidenceScores=ecModel.rxnConfidenceScores(index);
end
if isfield(ecModel,'pwys')
    rxnsToAdd.pwys=ecModel.pwys(index);
end

end