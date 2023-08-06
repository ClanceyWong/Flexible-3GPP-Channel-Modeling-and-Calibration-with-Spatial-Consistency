function SOS_model = sos_applygradient(SOS_model) % 都是各种优化方法，我们选normal

if strcmp(SOS_model.optimization_method, 'Adam')
    grad_squared = 0;
    for k = 1 : SOS_model.depth-1
        grad_squared = grad_squared + sum(sum(SOS_model.W_grad{k}.^2));
    end
end

for k = 1 : SOS_model.depth-1
    if strcmp(SOS_model.optimization_method, 'normal')
        if k~= SOS_model.depth-1
        SOS_model.W{k} = SOS_model.W{k} - SOS_model.learning_rate*SOS_model.W_grad{k};
        else
            if all((SOS_model.W{k} - SOS_model.learning_rate*SOS_model.W_grad{k})>0)
                SOS_model.W{k} = SOS_model.W{k} - SOS_model.learning_rate*SOS_model.W_grad{k};
                SOS_model.W{k} = SOS_model.W{k}/sum(SOS_model.W{k});
            end
        end
    elseif strcmp(SOS_model.optimization_method, 'Adam')
        rho1 = 0.9;
        rho2 = 0.999;
        delta = 0.00001;
        SOS_model.sW{k} = rho1*SOS_model.sW{k} + (1-rho1)*SOS_model.W_grad{k};
        SOS_model.rW{k} = rho2*SOS_model.rW{k} + (1-rho2)*SOS_model.W_grad{k}.^2;

        newS = SOS_model.sW{k}/(1-rho1^SOS_model.AdamTime);
        newR = SOS_model.rW{k}/(1-rho2^SOS_model.AdamTime);
        if k~= SOS_model.depth-1
            SOS_model.W{k} = SOS_model.W{k} - SOS_model.learning_rate*newS./sqrt(newR+delta);
        else
            if all((SOS_model.W{k} - SOS_model.learning_rate*newS./sqrt(newR+delta))>0)
                SOS_model.W{k} = SOS_model.W{k} - SOS_model.learning_rate*newS./sqrt(newR+delta);
                SOS_model.W{k} = SOS_model.W{k}/sum(SOS_model.W{k});
            end
        end
    end
end
end