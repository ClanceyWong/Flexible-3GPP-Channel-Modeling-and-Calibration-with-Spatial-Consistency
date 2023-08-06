function [error,SOS_model] = sos_test(SOS_model,test_xyz,test_acf)
    SOS_model = sos_predict(SOS_model,test_xyz);
    acf_output = SOS_model.a{SOS_model.depth};
    error = (1/numel(test_acf))*sum((acf_output' - test_acf).^2);
end