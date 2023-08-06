classdef blockage
    %BLOCKAGE 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        lambda
        sos
        SelfBlocking_LCS
        mode
        K = 0
        r
        nonSelfBlocking_GCS
        scenario
        v = 0
    end
    
    methods
        function obj = blockage(varargin)
            if isempty(varargin)
                obj.lambda   = 3e8/3e10;
                obj.mode     = 'Landscape';
                obj.scenario = '3D-UMi';
            else
                obj.lambda   = 3e8/varargin{1};
                obj.mode     = varargin{2};
                obj.scenario = varargin{3};
            end
            if strcmp(obj.scenario,'3D-InH')
                obj.sos{1}      = tools.sos('exp_300');
                obj.sos{2}      = tools.sos('exp_300');
                d_corr = [5,5];
                obj.sos{1}.init(repmat(d_corr(1),1,obj.K*3));
                obj.sos{2}.init(repmat(d_corr(2),1,obj.K*3));
            elseif strcmp(obj.scenario,'3D-UMi') ||strcmp(obj.scenario,'3D-UMa')||strcmp(obj.scenario,'3D-RMa')
                obj.sos{1}      = tools.sos('exp_300');
                obj.sos{2}      = tools.sos('exp_300');
                obj.sos{3}      = tools.sos('exp_300');
                d_corr = [10,10,5];
                obj.sos{1}.init(repmat(d_corr(1),1,obj.K*2));
                obj.sos{2}.init(repmat(d_corr(2),1,obj.K*2));
                obj.sos{3}.init(repmat(d_corr(3),1,obj.K*2));
            else
                error("Only support the scenario of '3D-InH', '3D-UMi', '3D-UMa' or '3D-RMa'.");
            end
%             if obj.v ~= 0
%                 obj.sos{2}      = tools.sos('exp_300');
%                 obj.sos{2}.init(d_corr./v);
%             end
            obj = obj.getSelfBlocking_region_parameters();
        end
        
        function [lossSelf, loss] = attenuation(obj,thetaLCS,phiLCS,thetaGCS,phiGCS,pos,indx)
            lossSelf = 0; loss = 0;
            phiLCS(phiLCS<0) = phiLCS + 360;
            phiGCS(phiGCS<0) = phiGCS + 360;
            if abs(thetaLCS-obj.SelfBlocking_LCS(3))<(obj.SelfBlocking_LCS(4)/2) && min(abs(phiLCS-obj.SelfBlocking_LCS(1)),360-abs(phiLCS-obj.SelfBlocking_LCS(1)))<(obj.SelfBlocking_LCS(2)/2)
                lossSelf = -30;% dB
            end
            obj = obj.getBlocking_region_parameters(pos,indx);
            for k = 1:obj.K
                if abs(thetaGCS-obj.nonSelfBlocking_GCS(k,3))<(obj.nonSelfBlocking_GCS(k,4)) && min(abs(phiGCS-obj.nonSelfBlocking_GCS(k,1)),360-abs(phiGCS-obj.nonSelfBlocking_GCS(1)))<(obj.nonSelfBlocking_GCS(k,2))
                    F = obj.F_AZ12(thetaGCS,phiGCS,obj.nonSelfBlocking_GCS(k,:));
                    loss = 20*log10(1-(F(1)+F(2))*(F(3)+F(4))) + loss;
                end
            end
        end
        
        function F = F_AZ12(obj,thetaGCS,phiGCS,nonSelfBlocking_GCS)
            sign = [1,1,1,1];
            if thetaGCS-nonSelfBlocking_GCS(3)<=-nonSelfBlocking_GCS(4)/2
                sign(4) = -1;
            elseif thetaGCS-nonSelfBlocking_GCS(3)>nonSelfBlocking_GCS(4)/2
                sign(3) = -1;
            end
            if phiGCS-nonSelfBlocking_GCS(1)>0 || phiGCS-nonSelfBlocking_GCS(1)<-180
                phi = min(abs(phiGCS-nonSelfBlocking_GCS(1)),360-abs(phiGCS-nonSelfBlocking_GCS(1)));
            elseif phiGCS-nonSelfBlocking_GCS(1)<0 || phiGCS-nonSelfBlocking_GCS(1)>180
                phi = max(-(phiGCS-nonSelfBlocking_GCS(1)),-360+abs(phiGCS-nonSelfBlocking_GCS(1)));
            end
            if phi>nonSelfBlocking_GCS(2)/2 
                sign(1) = -1;
            elseif phi<=-nonSelfBlocking_GCS(2)/2
                sign(2) = -1;
            end
            AZ12 = [phi-nonSelfBlocking_GCS(2)/2,...
                    phi+nonSelfBlocking_GCS(2)/2,...
                    thetaGCS-nonSelfBlocking_GCS(3)-nonSelfBlocking_GCS(4)/2,...
                    thetaGCS-nonSelfBlocking_GCS(3)+nonSelfBlocking_GCS(4)/2];
            
            F = (atan(sign*pi/2.*sqrt(pi/obj.lambda*obj.r*((1./cosd(AZ12))-1))))/pi;
        end
        
        function obj = getSelfBlocking_region_parameters(obj,varargin) % mode = 'Portrait' or 'Landscape'
            if isempty(varargin)
                mode_ = obj.mode;
            else
                mode_ = varargin{1};
            end
            if strcmp(mode_,'Portrait')
                params = [260;120;100;80]; % phi,x_sb,theta,y_sb;
            elseif strcmp(mode_,'Landscape')
                params = [40;160;110;75];  % phi,x_sb,theta,y_sb;
            else
                error("Only support the mode of 'Portrait' or 'Landscape'.");
            end
            obj.SelfBlocking_LCS = params.';
        end
        
        function obj = getBlocking_region_parameters(obj,pos,indx)
            scenario_ = obj.scenario;
            randnum = obj.sos{indx}.rand(pos).';
            if strcmp(scenario_,'3D-InH')
                obj.r = 2;
                params = [(360*randnum(1:obj.K));(30*randnum(obj.K+1:2*obj.K)+15);90*ones(1,obj.K);(10*randnum(2*obj.K+1:3*obj.K)+5)]; % phi,x_sb,theta,y_sb;
            elseif strcmp(scenario_,'3D-UMi') ||strcmp(scenario_,'3D-UMa')||strcmp(scenario_,'3D-RMa')
                obj.r = 10;
                params = [(360*randnum(1:obj.K));(10*randnum(obj.K+1:2*obj.K)+5);90*ones(1,obj.K);5*ones(1,obj.K)];  % phi,x_sb,theta,y_sb;
            else
                error("Only support the scenario of '3D-InH', '3D-UMi', '3D-UMa' or '3D-RMa'.");
            end
            obj.nonSelfBlocking_GCS = params.';
        end
    end
end

