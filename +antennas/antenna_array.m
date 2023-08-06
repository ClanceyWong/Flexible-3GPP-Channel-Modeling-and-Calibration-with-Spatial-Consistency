classdef antenna_array < handle
    % ANTENNA_ARRAY
    
    properties
        Mg                     % Mg rows
        Ng                     % Ng columns
        dg_H                   % Horizontal antenna element spacing
        dg_V                   % Vertical antenna element spacing 
        
        alpha                 % bearing angle for LCS to GCS
        beta                  % downtilt angle for LCS to GCS, (boresight)
        gamma                 % slant angle for LCS to GCS
        
        attachedType          % 'BS' or 'UE'
        attachedDevice        % attached Sector or UE
        
        panel                 % [Mg x Ng] class {antenna_panel}
    end
    properties(Dependent)
        num_panel
        pos_panel_LCS
        R
    end
    
    methods
        function obj = antenna_array(varargin)
            % input: 
            %        varargin{1} : Mg, Ng, dg_H, dg_V;
            %        varargin{2} : alpha, beta, gamma
            %        varargin{3} : a class of {antennas.antenna_panel}
            % output:
            %        obj : the class of {antennas.antenna_array}
            if isempty(varargin) 
                % default parameters
                obj.Mg     = 1;
                obj.Ng     = 1;
                obj.dg_H   = 2.5;
                obj.dg_V   = 2.5;
                obj.alpha = 0; 
                obj.beta  = 0; 
                obj.gamma = 0;
                obj.panel = antennas.antenna_panel();
                obj.panel.attachedArray = obj;
                obj.panel.ID = [1,1];
            else
                obj.Mg     = varargin{1}.array.Mg; 
                obj.Ng     = varargin{1}.array.Ng; 
                obj.dg_H   = varargin{1}.array.dg_H; 
                obj.dg_V   = varargin{1}.array.dg_V; 
                
                obj.alpha = varargin{2}.alpha; 
                obj.beta  = varargin{2}.beta; 
                obj.gamma = varargin{2}.gamma;
                
                obj.panel = antennas.antenna_panel;
                for m = 1:obj.Mg
                    for n = 1:obj.Ng
                        obj.panel(m,n)    = antennas.antenna_panel(varargin{1}.panel);
                        obj.panel(m,n).attachedArray = obj;
                        obj.panel(m,n).ID = [m,n];
                    end
                end
            end
        end
        
        function gain = array_gain(obj,link)
            gain = cell(obj.Mg,obj.Ng);
            for m = 1:obj.Mg
                for n = 1:obj.Ng
                    gain{m,n} = obj.panel(m,n).panel_gain(link);
                end
            end
        end
        
        function set_num_panel(obj,mn) % mg,ng
            % set the number of antenna_panel
            if ~(all(size(mn) == [1,2]) || all(size(mn) == [2,1]))
                error('Please input the parameters with the size of [2 x 1] or [1 x 2]');
            end
            obj.Mg        = mn(1);
            obj.Ng        = mn(2);
        end
        
        %% get and set function
        function out = get.num_panel(obj)
            % get the number of antenna_panel
            out = obj.Mg*obj.Ng;
        end
        function out = get.pos_panel_LCS(obj)
            % get the position of the antenna panel in LCS
            [m, n]     = meshgrid(0:obj.Mg-1,0:obj.Ng-1);
            out        = zeros(obj.Mg,obj.Ng,3);
            out(:,:,2) = (n'-(obj.Ng-1)/2)*obj.dg_H;
            out(:,:,3) = (m'-(obj.Mg-1)/2)*obj.dg_V;
        end
        function out = get.R(obj)
            a = obj.alpha;
            b = obj.beta;
            c = obj.gamma;
            out  = [cosd(a)*cosd(b)                        , sind(a)*cosd(b)                        , -sind(b);
                    cosd(a)*sind(b)*sind(c)-sind(a)*cosd(c), sind(a)*sind(b)*sind(c)+cosd(a)*cosd(c), cosd(b)*sind(c);
                    cosd(a)*sind(b)*cosd(c)+sind(a)*sind(c), sind(a)*sind(b)*cosd(c)-cosd(a)*sind(c), cosd(b)*cosd(c)].';
            
        end
    end
end