function SOS_model = sos_create(varargin)%varargin�൱��c��argv[]
    %�ṹ�����ڸ��ֲ���������ѡ������
    SIZE = varargin{1};
    SOS_model.keep_probability    = 1;
    SOS_model.size                = SIZE;
    SOS_model.depth               = numel(SIZE); % ���������SIZEԪ�ظ�����
    SOS_model.active_function     = 'cos';
    SOS_model.output_function     = 'sum';
    SOS_model.learning_rate       = 1.5; % ѧϰ��
    SOS_model.weight_decay        = 0; % Ϊ��ֹ����ϣ�weight_decay�Ƿ���������ǰ���һ��ϵ��
    SOS_model.cost                = [];
    SOS_model.grad_squared        = 0;
    SOS_model.r                   = 0;
    SOS_model.optimization_method = 'normal';
    SOS_model.objective_function  = 'MSE'; % ����������Error
    SOS_model.ad{1}               = 0;

    for i = 2:length(varargin)
        if strcmp('ActiveFunction',varargin{i})
            SOS_model.active_function = varargin{i+1};
        elseif strcmp('OutputFunction',varargin{i})
            SOS_model.output_function = varargin{i+1};
        elseif strcmp('LearningRate',varargin{i})
            SOS_model.learning_rate = varargin{i+1};
        elseif strcmp('WeightDecay',varargin{i})
            SOS_model.weight_decay = varargin{i+1};
        elseif strcmp('Sparsity',varargin{i})
            SOS_model.sparsity = varargin{i+1};
        elseif strcmp('OptimizationMethod',varargin{i})
            SOS_model.optimization_method = varargin{i+1};
        elseif strcmp('ObjectiveFunction', varargin{i})
            SOS_model.objective_function = varargin{i+1};
        elseif strcmp('KeepProbability',varargin{i})
            SOS_model.keep_probability = varargin{i+1};
        end
    end

    for k = 1 : SOS_model.depth - 1
        width_in = SOS_model.size(k);
        height_out = SOS_model.size(k+1);

        if k == SOS_model.depth-1
            SOS_model.W{k} = ones(height_out, width_in)/width_in;%rand����α��������󣬼�WȨ�ؾ����ʼ��
        else
            [ theta, phi ] = fibonacci_sphere( height_out );
            fn = ( randn(height_out,1)*1);
            if width_in == 3
                SOS_model.W{k} = [ fn .* cos( phi ) .* cos( theta ),...
                    fn .* sin( phi ) .* cos( theta ),...
                    fn .* sin( theta )]; % ����α��������󣬼�WȨ�ؾ����ʼ��
            elseif width_in == 6
                SOS_model.W{k} = [ fn .* cos( phi ) .* cos( theta ),...
                    fn .* sin( phi ) .* cos( theta ),...
                    fn .* sin( theta ),...
                    fn .* cos( phi ) .* cos( theta ),...
                    fn .* sin( phi ) .* cos( theta ),...
                    fn .* sin( theta )];
            end
        end
        if abs(SOS_model.keep_probability-1)>0.001
            SOS_model.WMask{k} = ones(height_out,width_in);
        end

        % initialization
        SOS_model.b{k} = zeros(height_out, 1); % b��ֵ�ĳ�ʼ��

        % parameters for optimization method
        if strcmp(SOS_model.optimization_method,'Adam')
            SOS_model.rW{k} = zeros(height_out,width_in);
            SOS_model.sW{k} = zeros(height_out,width_in);
            % sos.rb{k} = zeros(height,1);
            % sos.sb{k} = zeros(height,1);
        end
        SOS_model.W_grad{k} = zeros(height_out,width_in);
    end
    if  strcmp(SOS_model.optimization_method,'Adam')
        SOS_model.AdamTime = 0;
    end
end
