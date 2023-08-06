classdef GroundReflection
    % GROUNDREFLECTION 此处显示有关此类的摘要
    %   此处显示详细说明
    properties(Constant)
        epsilon_0 = 8.854187817*1e-12;
        Concrete = 1
        Brick = 2
        Plasterboard = 3
        Wood = 4
        Floorboard = 5
        Metal = 6
        VeryDryGround = 7
        MediumDryGround = 8
        WetGround = 9
        list = [5.31	0	0.0326	0.8095;
                    3.75	0	0.038	0     ;
                    2.94	0	0.0116	0.7076;
                    1.99	0	0.0047	1.0718;
                    3.66	0	0.0044	1.3515;
                    1	    0	107	    0     ;
                    3	    0	0.00015	2.52  ;
                    15	 -0.1	0.035	1.63  ;
                    30	 -0.4	0.15	1.30];
    end
    properties
        Material
    end
    
    methods
        function obj = GroundReflection(material)
            obj.Material = material;
        end
        
        % Material properties
        function k = GetProperties(obj)
            k = obj.list(obj.Material,:);
        end
        
        % the complex relative permittivity of the ground material
        function out = GetRelaPerm(obj,fc)
            paras = obj.GetProperties;
            cond = paras(:,3)*(fc/10e9).^paras(:,4);
            perm = paras(:,1)*(fc/10e9).^paras(:,2);
            out  = perm - 1j*(cond./(2*pi*fc*obj.epsilon_0));
        end
        
        function Phi = GetPhi(obj,fc,theta_GR_ZOD)
            CRP = obj.GetRelaPerm(fc);
            Rpara = (CRP*cosd(theta_GR_ZOD) + sqrt(CRP-(sind(theta_GR_ZOD)).^2))/...
                    (CRP*cosd(theta_GR_ZOD) - sqrt(CRP-(sind(theta_GR_ZOD)).^2));
            Rperp = (cosd(theta_GR_ZOD) + sqrt(CRP-(sind(theta_GR_ZOD)).^2))/...
                    (cosd(theta_GR_ZOD) - sqrt(CRP-(sind(theta_GR_ZOD)).^2));
            Phi   = [Rpara, 0; 0, -Rperp];
        end
    end
end

