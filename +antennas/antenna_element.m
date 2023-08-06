classdef antenna_element < handle
    % ANTENNA_ELEMENT 
    
    properties
%         ID                    % antenna_element ID
        attachedPanel
%         pos_LCS               % relative 3D position
        slant_angle           % used for the polarized antennas, TR36.873 page 19
%         w_m          = 1      % Complex weight for antenna element m in elevation
        pol_model = 'model-2'
    end
    
    properties(Dependent)
        max_gain              % dBi
        attachedType          % 'BS' or 'UE'
        alpha                 % bearing angle for LCS to GCS
        beta                  % downtilt angle for LCS to GCS
        gamma                 % slant angle for LCS to GCS
        % POL                   % 'C' for vertically polarized, 'X' for cross-polarized
    end
    
    methods
        function obj = antenna_element(varargin)
            % antenna_element(slant_angle,alpha,beta,gamma)
            if isempty(varargin)
                obj.slant_angle = 0;
            else
                obj.slant_angle = varargin{1};
            end
%             obj.ID = ID;
%             obj.attachedType = attachedType;
%             obj.pos = pos;
        end
        
        function lcsgain = lcs_gain(obj, phi, theta) % LCS gain, phi and theta in LCS
            if strcmp(obj(1).attachedType,'BS') || isempty(obj(1).attachedType)
                Am       = 30; 
                phi3db   = 65;
                A_E_H    = -min(12*(phi/phi3db).^2, Am);
                SLAv     = 30; 
                theta3db = 65;
                A_E_V    = -min(12*((theta - 90)/theta3db).^2, SLAv);
                A        = -min(-(A_E_H+A_E_V), Am);
                lcsgain  = 10.^((A+obj.max_gain)/10);  % in power
            elseif strcmp(obj(1).attachedType,'UE')
                lcsgain  = ones(size(phi)); % in power
            elseif strcmp(obj(1).attachedType,'BSInF')
                lcsgain  = ones(size(phi)); % in power
            elseif strcmp(obj(1).attachedType,'BS3D-InH')
                lcsgain  = ones(size(phi)); % in power
            end
        end
        
        function gcsgain = gcs_gain(obj, phi, theta) % GCS gain, phi and theta in GCS
            alpha_ = obj.alpha;
            beta_  = obj.beta;
            gamma_ = obj.gamma;
            inv_R  = [cosd(alpha_)*cosd(beta_)                                       , sind(alpha_)*cosd(beta_)                                       , -sind(beta_);
                      cosd(alpha_)*sind(beta_)*sind(gamma_)-sind(alpha_)*cosd(gamma_), sind(alpha_)*sind(beta_)*sind(gamma_)+cosd(alpha_)*cosd(gamma_), cosd(beta_)*sind(gamma_);
                      cosd(alpha_)*sind(beta_)*cosd(gamma_)+sind(alpha_)*sind(gamma_), sind(alpha_)*sind(beta_)*cosd(gamma_)-cosd(alpha_)*sind(gamma_), cosd(beta_)*cosd(gamma_)];
            xyz    = sph2cart_(phi, theta);
            tmp    = multiprod(inv_R,xyz,[1 2],1);
            tmp1   = tmp(1,:,:); 
            tmp2   = tmp(2,:,:); 
            tmp3   = tmp(3,:,:);
            [phi_tmp, theta_tmp] = cart2sph_(tmp1, tmp2, tmp3);
            gcsgain              = obj.lcs_gain(phi_tmp, theta_tmp);
        end
        
        function F_theta_phi = field_pattern(obj, phi, theta)
%             if isempty(varargin)
%                 pol_model = 'model-2';
%             else 
%                 pol_model = varargin{1};
%             end
            alpha_ = obj.alpha; beta_ = obj.beta; gamma_ = obj.gamma;
            if strcmp(obj.pol_model,'model-2')
                A_tmp          = obj.gcs_gain(phi, theta);
                F_theta_tmp    = sqrt(A_tmp)*cosd(obj.slant_angle);
                F_phi_tmp      = sqrt(A_tmp)*sind(obj.slant_angle);
                F_tmp(1,:,:,:) = F_theta_tmp; 
                F_tmp(2,:,:,:) = F_phi_tmp;
            elseif strcmp(obj.pol_model,'model-1')
                zeta            = obj.slant_angle;
                
                A_tmp2          = obj.gcs_gain(phi, theta);
                % antenna array rotation
                inv_R = [cosd(alpha_)*cosd(beta_)                                       , sind(alpha_)*cosd(beta_)                                       , -sind(beta_);
                         cosd(alpha_)*sind(beta_)*sind(gamma_)-sind(alpha_)*cosd(gamma_), sind(alpha_)*sind(beta_)*sind(gamma_)+cosd(alpha_)*cosd(gamma_), cosd(beta_)*sind(gamma_);
                         cosd(alpha_)*sind(beta_)*cosd(gamma_)+sind(alpha_)*sind(gamma_), sind(alpha_)*sind(beta_)*cosd(gamma_)-cosd(alpha_)*sind(gamma_), cosd(beta_)*cosd(gamma_)];

                xyz   = sph2cart_(phi, theta);
                tmp   = multiprod(inv_R',xyz,[1 2],1);
                tmp1  = tmp(1,:,:); tmp2 = tmp(2,:,:); tmp3 = tmp(3,:,:);
                [phi_tmp2, theta_tmp2] = cart2sph_(tmp1, tmp2, tmp3);

                F_theta_tmp2    = sqrt(A_tmp2);
                F_phi_tmp2      = zeros(size(F_theta_tmp2));
                cos_psi         = (cosd(zeta).*sind(theta_tmp2)+sind(zeta).*sind(phi_tmp2).*cosd(theta_tmp2))./sqrt(...
                                   1-(cosd(zeta).*cosd(theta_tmp2)-sind(zeta).*sind(phi_tmp2).*sind(theta_tmp2)).^2);
                sin_psi         = sind(zeta).*cosd(phi_tmp2)./sqrt(...
                                   1-(cosd(zeta).*cosd(theta_tmp2)-sind(zeta).*sind(phi_tmp2).*sind(theta_tmp2)).^2);
                rota(1,1,:,:)   = cos_psi; 
                rota(1,2,:,:)   = -sin_psi;
                rota(2,1,:,:)   = sin_psi; 
                rota(2,2,:,:)   =  cos_psi;
                F_tmp2(1,:,:,:) = F_theta_tmp2; 
                F_tmp2(2,:,:,:) = F_phi_tmp2;
                F_tmp           = multiprod(rota,F_tmp2,[1 2],[1 2]);
            else 
                error('only support model-1 and model-2');
            end
            cos_Psi = (cosd(beta_).*cosd(gamma_).*sind(theta)-(...
                        sind(beta_).*cosd(gamma_).*cosd(phi-alpha_)-sind(gamma_).*sind(phi-alpha_)).*cosd(theta))./sqrt(1-(...
                        cosd(beta_).*cosd(gamma_).*cosd(theta)+(...
                        sind(beta_).*cosd(gamma_).*cosd(phi-alpha_)-sind(gamma_).*sind(phi-alpha_)).*sind(theta)).^2);
            sin_Psi = (sind(beta_).*cosd(gamma_).*sind(phi-alpha_)+sind(gamma_).*cosd(phi-alpha_))./sqrt(1-(...
                        cosd(beta_).*cosd(gamma_).*cosd(theta)+(...
                        sind(beta_).*cosd(gamma_).*cosd(phi-alpha_)-sind(gamma_).*sind(phi-alpha_)).*sind(theta)).^2);
            rota(1,1,:,:) = cos_Psi;  rota(1,2,:,:) = -sin_Psi;
            rota(2,1,:,:) = sin_Psi;  rota(2,2,:,:) =  cos_Psi;
            F_theta_phi   = multiprod(rota,F_tmp,[1 2],[1 2]);
        end
        
        function plt(obj) % plot the antenna element pattern
            figure;
            subplot(221);
            g1 = 10*log10(obj.lcs_gain(0,(0:1:180)));
            plot((0:1:180),g1,'linewIDth',1.5);
            grid on; xlim([0 180]);
            xlabel('Zenith angle \theta (^o)');ylabel('Gain (dB)');
            subplot(222);
            g2 = 10*log10(obj.lcs_gain((-180:1:180),90));
            plot((-180:1:180),g2,'linewIDth',1.5);
            grid on; xlim([-180 180]);
            xlabel('Azimuth angle \phi (^o)');ylabel('Gain (dB)');
            
            theta         = linspace(0,pi,200);
            phi           = linspace(0,2*pi-0.0001,200);
            [phi_,theta_] = meshgrid(phi,theta);
            r             = obj.lcs_gain((phi_*180/pi-360*floor(phi_/pi)),(theta_*180/pi)); % in power
            [x,y,z]       = sph2cart((phi_*180/pi-360*floor((phi_-eps)/pi))*pi/180,(pi/2-theta_),r);
            subplot(2,2,[3 4]);
            mesh(x,y,z);
            title('Pattern of TR36.873 3D Antenna element');
        end
        
        %% get and set function
        function out = get.max_gain(obj)
            if strcmp(obj.attachedType(1:2),'BS')||isempty(obj.attachedType)
                out = 8;
            elseif strcmp(obj.attachedType,'UE')
                out = 1;
            else
                error('Unknown attached type of device.');
            end
        end
        function out = get.attachedType(obj)
            out = obj.attachedPanel.attachedType;            
        end
        function out = get.alpha(obj)
            out = obj.attachedPanel.alpha;            
        end
        function out = get.beta(obj)
            out = obj.attachedPanel.beta;            
        end
        function out = get.gamma(obj)
            out = obj.attachedPanel.gamma;            
        end
    end
end

% other function
function [phi, theta] = cart2sph_(x, y, z)
    [phi,ele,~] = cart2sph(x,y,z);
    phi         = phi*180/pi;          % [-180, 180]
    theta       = 90 - ele*180/pi;     % [   0, 180]
end
function [xyz] = sph2cart_(phi, theta, varargin)
    if isempty(varargin)
        r = ones(size(phi));
    else 
        r = varargin{1};
    end
    x              = r.*sind(theta).*cosd(phi);
    y              = r.*sind(theta).*sind(phi);
    z              = r.*cosd(theta);
    [n, m]         = size(phi);
    xyz(1,1:n,1:m) = x;
    xyz(2,1:n,1:m) = y;
    xyz(3,1:n,1:m) = z;
%     xyz = [x y z];
end