function SOS_model = sos_predict(SOS_model,batch_xyz) % 前面文件中都注释过
    batch_xyz = batch_xyz';
    m = size(batch_xyz,2);
    SOS_model.a{1} = batch_xyz;
    for k = 2 : SOS_model.depth
        y = SOS_model.W{k-1} * SOS_model.a{k-1} + repmat(SOS_model.b{k-1},1,m);
        if k == SOS_model.depth
            switch SOS_model.output_function
                case 'sigmoid'
                    SOS_model.a{k} = sigmoid(y);
                case 'tanh'
                    SOS_model.a{k} = tanh(y);
                case 'relu'
                    SOS_model.a{k} = max(y,0);
                case 'softmax'
                    SOS_model.a{k} = softmax(y);
                case 'sum'
                    SOS_model.a{k} = y;
            end
        else
            switch SOS_model.active_function
                case 'sigmoid'
                    SOS_model.a{k} = sigmoid(y);
                case 'tanh'
                    SOS_model.a{k} = tanh(y);
                case 'relu'
                    SOS_model.a{k} = max(y,0);
                case 'cos'
                    SOS_model.a{k} = cos(y);
            end
        end
    end
end