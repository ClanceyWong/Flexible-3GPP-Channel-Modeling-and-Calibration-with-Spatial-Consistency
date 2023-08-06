function SOS_model = sos_train(SOS_model,option,train_xyz,train_acf)
    iteration = option.iteration;     % ��������������������
    batch_size = option.batch_size;   % �ۻ�BP�㷨�е�������С
    m = size(train_xyz,1);            % size(��1)���ؾ�������
    num_batches = m / batch_size;
    for k = 1 : iteration
        index = randperm(m);             % ����1��m��˳�� 
        if  strcmp(SOS_model.optimization_method,'Adam')
            SOS_model.AdamTime = SOS_model.AdamTime+1;
        end
        for batche = 1 : num_batches       % ����60000��ѵ�����ݴ��Һ���С��
            batch_xyz = train_xyz(index((batche-1)*batch_size + (1:batch_size)),:);
            batch_acf = train_acf(index((batche-1)*batch_size + (1:batch_size)),:);
            SOS_model = sos_forward(SOS_model,batch_xyz,batch_acf); % ǰ�����
            SOS_model = sos_backpropagation(SOS_model,batch_acf);   % �������
            SOS_model = sos_applygradient(SOS_model);               % �ݶ��½�
        end
    end
end