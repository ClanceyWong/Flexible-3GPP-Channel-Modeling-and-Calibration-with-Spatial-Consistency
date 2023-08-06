function SOS_model = sos_backpropagation(SOS_model,batch_acf)
    batch_acf = batch_acf';
    m = size(SOS_model.a{1},2);  % m Ϊ batch ����Ŀ
    SOS_model.theta{1} = 0;      % ������¼ÿһ��ĵ�����С����Ҳ�Ǻ��򴫲�������
    switch SOS_model.output_function
        case 'sigmoid'
            SOS_model.theta{SOS_model.depth} = -(batch_acf-SOS_model.a{SOS_model.depth}) .* SOS_model.a{SOS_model.depth} .* (1 - SOS_model.a{SOS_model.depth});
        case 'tanh'
            SOS_model.theta{SOS_model.depth} = -(batch_acf-SOS_model.a{SOS_model.depth}) .* (1 - SOS_model.a{SOS_model.depth}.^2);
        case 'softmax'
            SOS_model.theta{SOS_model.depth} = SOS_model.a{SOS_model.depth} - batch_acf;
        case 'sum'
            SOS_model.theta{SOS_model.depth} = SOS_model.a{SOS_model.depth} - batch_acf;
    end

    % �ڶ�����������������󵼲����ģ���������ϵ�� SOS_model.weight_decay ���� w
    SOS_model.W_grad{SOS_model.depth-1} = SOS_model.theta{SOS_model.depth}*SOS_model.a{SOS_model.depth-1}'/m + SOS_model.weight_decay*SOS_model.W{SOS_model.depth-1};
    % ����ʽ�Ӷ��� m ���ݶ�����ƽ��
    switch SOS_model.active_function
        case 'cos'
            for ll = 2 : SOS_model.depth - 1
                k = SOS_model.depth - ll + 1; % ������õ��˺��򴫲��ĵ������壬�Ӻ���ǰ����
                SOS_model.theta{k} = ((SOS_model.W{k}'*SOS_model.theta{k+1})) .* SOS_model.ad{k};
                SOS_model.W_grad{k-1} = SOS_model.theta{k}*SOS_model.a{k-1}'/m + SOS_model.weight_decay*SOS_model.W{k-1};
            end
    end
end