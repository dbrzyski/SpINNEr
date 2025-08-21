classdef spinner < handle
    properties
        y
        X
        A % matrix of a form [A1,A2,...]
        AA % 3- way tensor
        W

        % sizes
        n
        p
        d
        % Params
        % ParamsFit
        % usedDefaults
        % usedDefaultsFit
    end

    methods
        %------------------------------------------------------------------
        % Initialization
        %------------------------------------------------------------------
        function obj = spinner(y, X, A, W)
            if nargin == 0
                % show all the list of all possible arguments ( and the corresponding default values)
                defaultSettings = ParseArgumentsFit({});
                disp(defaultSettings)
                return;
            elseif nargin < 4
                StandardWeightsMatrix = true; % 'Zeros on diagonal, ones on off-diagonal'
            else
                StandardWeightsMatrix = false;
            end

            % Assign the required arguments
            obj.y = y;
            obj.X = X;
            obj.d = size(X,2);

            % Tensor of regressor matrices
            obj.convertRegressorMatrices(A);

            % Collect input arguments
            % [obj.Params, obj.usedDefaults] = ParseArguments(varargin);

            if StandardWeightsMatrix
                obj.W = ones(obj.p, obj.p) - eye(obj.p);
            else
                obj.W = W;
            end
        end

        %------------------------------------------------------------------
        function convertRegressorMatrices(obj, A)
            [p1, p2, p3] = size(A);
            if p3==1
                if mod(p2,p1)~=0
                    error('The assumed form is AA = [A1, A2,...,An], with square matrices Ai, hence number of columns in AA should be the multiplicity of the number of rows')
                else
                    obj.A = A;
                    obj.AA = reshape(A, [p1, p1, p2/p1]);
                end
            else
                if p1~=p2
                    error('The matrices Ai are assumed to be square')
                else
                    obj.A = reshape(A, [p1*p2, p3])';
                    obj.AA = A;
                end
            end

            % sizes
            obj.p = p1;
            obj.n = p3;
        end

        %------------------------------------------------------------------
        function out = fit(obj, varargin)
            % if nargin == 1
            %     % show all the list of all possible arguments ( and the corresponding default values)
            %     defaultSettings = ParseArgumentsFit(varargin);
            %     disp(defaultSettings)
            %     return;
            % end
            % Collect input arguments
            [ParamsFit, ~] = ParseArgumentsFit(varargin);

            % % extract the information on W and remove it from ParamFit
            % if isequal(ParamsFit.W, 'Zeros on diagonal, ones on off-diagonal')
            %     cW = ones(obj.p, obj.p) - eye(obj.p);
            % else
            %     cW = ParamsFit.W;
            % end
            % ParamsFit = rmfield(ParamsFit, 'W');

            %--------------------------------------------------------------
            %                  PROVIDED TUNING PARAMETERS
            %--------------------------------------------------------------

            if ~isempty(ParamsFit.LambdaL) || ~isempty(ParamsFit.LambdaN)
                if isempty(ParamsFit.LambdaL) || isempty(ParamsFit.LambdaN)
                    error('[not yet implemented] both or none parameters must be provided')
                else
                    % two parameters provided
                    if isequal(ParamsFit.Family, 'Gaussian')
                        out = spinnerRun(obj.y, obj.X, obj.AA, ParamsFit.LambdaN, ParamsFit.LambdaL, obj.W, ParamsFit);
                    elseif isequal(ParamsFit.Family, 'Binomial')
                        assert(all(or(obj.y==0, obj.y == 1)), 'Response vector, y, should have only zeros and ones coefficients if Family == "Binomial"')
                        out = LogisticSpinner(obj.y, obj.X, obj.AA, ParamsFit.LambdaN, ParamsFit.LambdaL, obj.W, ParamsFit);
                    else
                        error('Not implemented')
                    end
                end
            else

            %--------------------------------------------------------------
            %                 NO TUNING PARAMETERS
            %--------------------------------------------------------------

                if isequal(ParamsFit.Family, 'Gaussian')
                    if isequal(ParamsFit.Method, 'CV')
                        % % SPINNER: Find parameters by Cross-Validation
                        % cParams = {'UseParallel', 'gridLengthN', 'gridLengthL', 'gridParameter', ...
                        %     'kfolds', 'displayStatus', 'initLambda', 'zeroSearchRatio', 'maxLambAcc'};
                        % cArgs = obj.getSubstruct(ParamsFit, cParams);
                        % cArgs.W = cW;
                        out = spinnerCV(obj.y, obj.X, obj.AA, obj.W, ParamsFit);
                    else
                        % SPINNER: Find parameters by Bayesian
                        error('Not implemented')
                    end
                elseif isequal(ParamsFit.Family, 'Binomial')
                    assert(all(or(obj.y==0, obj.y == 1)), 'Response vector, y, should have only zeros and ones coefficients if Family == "Binomial"')
                    if isequal(ParamsFit.Method, 'CV')
                        % LOGISTIC SPINNER: Find parameters by Cross-Validation
                        % cParams = {'UseParallel', 'gridLengthN', 'gridLengthL', 'gridParameter',...
                        %     'kfolds', 'displayStatus', 'initLambda', 'zeroSearchRatio', 'maxLambAcc'};
                        % cArgs = obj.getSubstruct(ParamsFit, cParams);
                        % cArgs.W = cW;
                        out = Logistic_spinnerCV(obj.y, obj.X, obj.AA, obj.W, ParamsFit);
                    else
                        error('Not implemented')
                    end
                else
                    error('Not implemented')
                end
            end
        end
        
        %------------------------------------------------------------------
        function updateWeights(obj, W)
            % obj.OptionalArgs.W = W;
            obj.W = W;
        end

        %------------------------------------------------------------------
        function subStruct = getSubstruct(~, S, fields)
            Fields = fieldnames(S);
            foundFieldsIdxs = ismember(fields,Fields);
            if ~all(foundFieldsIdxs)
                missingFields = join(fields(~foundFieldsIdxs),', ');
                error(['The following fields names were not identified in struct: ', missingFields{1}])
            end
            
            nFields = length(fields);
            for ii = 1:nFields
                if ismember(fields{ii}, Fields)
                    subStruct.(fields{ii}) = S.(fields{ii});
                end
            end
        end
    end
end
