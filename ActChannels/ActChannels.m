classdef ActChannels < handle
    % class ActChannels to test basic kinetics of active channels. 
    % channels density is not considered here. 
    % todo: Ca dependent channels
    properties (SetAccess = protected)
        Vm;
        dt;
        ds;
        Ca1;
        Ca2;
        ICa;
        nch;
        nf;
        ECa;
        fchan;
        channelstype;
        ch_inf;
        ch_tau;
        channame;
    end
    methods
        function self=ActChannels(channelstype, dt, Vm)
            % initiation of channels 
            % self.channelstype=channelstype;
            % example: 
            %   actCh=ActChannels(channelT, [], []);
            % initiation will give the properties of all nominated channels
            self.dt = dt;
            if isempty(Vm)
                self.Vm = linspace(-90, 20, 200)';
            else
                self.Vm = Vm;
            end
           self.nch = length(channelstype);
            channames={'sdNa';...% 1 Somato-dendritic Na+ channel (Na_7)
                'aNa';...% 2 Axonic Na+ channel (NaA_2)
                'wCa';...% 3 High thereshold Ca++ channel (Ca_W)
                'tCa';...% 4 Low thereshold Ca++ channel (CaT_3)
                'aKp';...% 5 A type K+ channel, proximal (K_A_11)
                'aKd';...% 6 A type K+ channel, distal (K_A_18)
                'drK';...% 7 DR type K+ channel (K_DR_2)
                % 'drKa';...% DR type K+ channel, axonic (K_DRA_4)
                'cK';...% 8 C type K+ channel (K_C_1) ca 
                'mK';...% 9 M type K+ channel (K_M_4)
                'AHPK';...% 10 AHP type K+ channel (K_AHP_Wtn) ca
                'Ih';...%  h type channel (Ih_3)
                'UDefine';...% User defined channels
                };
            self.channame = cell(self.nch,1);
            self.fchan = cell(self.nch,1);
            self.ch_inf = cell(self.nch,2);% openning/closing
            self.ch_tau = cell(self.nch,2);% openning/closing
            
            for k=1:self.nch
                channame=channames{channelstype(k)};
                self.fchan{k}=str2func(channame);
                self=self.fchan{k}(self,k);
                self.channame{k}=channame;
            end
            self.plotAC
        end

        %% ---------------Somato-dendritic Na+ channel (Na_7)------- 1
        function self=sdNa(self,k)
                minf = 1./(1+exp((25+self.Vm)/-5));
                taum = .020+.22./(exp((self.Vm+55)/13)+exp((self.Vm+55)/-80));% ms
                self.ch_inf{k,1} = minf;
                self.ch_tau{k,1} = taum;
                
                hinf = 1./(1+exp((50+self.Vm)/4));
                tauh = .6+8*exp((self.Vm+45)/13)./(1+exp((self.Vm+45)/5.5));
                self.ch_inf{k,2}=hinf;
                self.ch_tau{k,2}=tauh;
        end
        %% ---------------------Axonic Na+ channel (NaA_2)------- 2
        function self=aNa(self,k)
                minf = 1./(1+exp((32+self.Vm)/-4.2));
                taum = .020+.220./(exp((self.Vm+62)/13)+exp((self.Vm+62)/-80));
                self.ch_inf{k,1} = minf;
                self.ch_tau{k,1} = taum;
                
                hinf = 1./(1+exp((50+self.Vm)/4));
                tauh = .0200+8*exp((self.Vm+45)/13)./(1+exp((self.Vm+45)/5.5));
                self.ch_inf{k,2}=hinf;
                self.ch_tau{k,2}=tauh;
        end
        
        %% -------------High thereshold Ca++ channel (Ca_W)----------- 3
        function self=wCa(self,k)
            self.Vm(abs(self.Vm-26)<10^-4 |abs(self.Vm-12)<10^-4 )=self.Vm(abs(self.Vm-26)<10^-4 |abs(self.Vm-12)<10^-4 )+10^-3;
            alpha = .160*(self.Vm+26)./(-exp((self.Vm+26)/-4.5)+1);
            beta = .040*(self.Vm+12)./(exp((self.Vm+12)/10)-1);
            minf=alpha./(alpha-beta);
            self.ch_inf{k,1} = minf;
            self.ch_tau{k,1} = 1./(alpha-beta);
            
            alpha = 2./exp((self.Vm+94)/10);
            beta = 8./(exp((self.Vm-68)/-27)+1);
            hinf=alpha./(alpha-beta);
            self.ch_inf{k,2}=hinf;
            self.ch_tau{k,2}=1./(alpha-beta);
            
        end
        %% -------------------Low thereshold Ca++ channel (CaT_3)---- 4
        function self=tCa(self,k)
            sinf = 1./(1+exp((54+self.Vm)/-5));
            taus = .204+.333./(exp((self.Vm+15.8)/18.2)+exp((self.Vm+131)/-16.7));
            self.ch_inf{k,1} = sinf;
            self.ch_tau{k,1} = taus;
            
            rinf = 1./(1+exp((80+self.Vm)/4));
            taur=rinf;
            taur(self.Vm<= -81)=.333*exp((self.Vm(self.Vm<= -81)+466)/66.6);
            taur(self.Vm>-81)=9.32 + .333*exp((self.Vm(self.Vm> -81)+21)/-10.5);
            self.ch_inf{k,2} = rinf;
            self.ch_tau{k,2} = taur;
        end
        %% --------------A type K+ channel, proximal (K_A_11)--------- 5
        function self=aKp(self,k)
            ninf = 1./(1+exp((12+self.Vm)/-8.5));
            taun = repmat(.080, size(self.Vm));
            self.ch_inf{k,1} = ninf;
            self.ch_tau{k,1} = taun;
            
            linf = 1./(1+exp((56+self.Vm)/7));
            taul = max(.26*(self.Vm+50),2);
            self.ch_inf{k,2} = linf;
            self.ch_tau{k,2} = taul;
        end
        
        %% ------------- A type K+ channel, distal (K_A_18)------ 6
        function self=aKd(self,k)
            ninf = 10.^(self.Vm<= -65)./(1+exp((24+self.Vm)/-7));
            taun =repmat(.080, size(self.Vm)); 
            self.ch_inf{k,1} = ninf;
            self.ch_tau{k,1} = taun;
            
            linf = 1./(1+exp((56+self.Vm)/7));
            taul = max(0.26*(self.Vm+50),2);
            self.ch_inf{k,2} = linf;
            self.ch_tau{k,2} = taul;
        end
        
        %% -------------DR type K+ channel (K_DR_2)--------------- 7
        function self=drK(self,k)
            ninf = 1./(1+exp((5-self.Vm)/11));
            taun = repmat(1.2, size(self.Vm));
            self.ch_inf{k,1} = ninf;
            self.ch_tau{k,1} = taun;
        end
        
        %% -------------C type K+ channel (K_C_1)-------------- 8
        function self=cK(self,k)
            % Ca dependent C type K+ channel
            Vshift = 40.*log(max((self.Ca1+.05)*1000,10^-7))-105;
            alpha = -7.7*(self.Vm+Vshift+103)/1000./(exp((self.Vm+Vshift+103)/-12)-1);
            beta = 1.7./exp((self.Vm+Vshift+237)/30);
            minf=alpha./(alpha-beta);
            self.ch_inf{k,1} = minf;
            self.ch_tau{k,1} = 1./(alpha-beta);
            
            
            alpha = 1./exp((self.Vm+79)/10);
            beta =4./(exp((self.Vm-82)/-27)+1);
            hinf=alpha./(alpha-beta);
            self.ch_inf{k,2} = hinf;
            self.ch_tau{k,2} = 1./(alpha-beta);
        end
        
        %% --------------- M type K+ channel (K_M_4)--------------- 9
        function self=mK(self,k)
            ninf = 1./(1+exp((40+self.Vm)/-8.5));
            taun =0.9+20*exp((self.Vm+38)/25)./(1+exp((self.Vm+38)/7));
            self.ch_inf{k,1} = ninf;
            self.ch_tau{k,1} = taun;
        end
        %% -------------AHP type K+ channel (K_AHP_Wtn)------------- 10
        function self=AHPK(self,k)
            % self=AHPK(self,k): Ca dependent AHP type K+ channel
            alpha = 4.8*10^6./exp(-5*log(max((self.Ca2+.05)*1000-35,10^-6)));
            beta = 12*10^6./exp(2*log(max((self.Ca2+.05)*1000+100,10^-6)));
            minf=alpha./(alpha-beta);
            self.ch_inf{k,1} = minf;
            self.ch_tau{k,1} = 1./(alpha-beta);
        end
        function h=plotAC(self)
            h=figure;
            for k=1:self.nch
                subplot(2, self.nch, k)
                plot(self.Vm, cell2mat(self.ch_inf(k,:)))
                hold on;
                plot([-70, -70], [0 1],'r--')
                grid on; axis tight
                try
                legend('C-O','O-C')
                catch
                legend('C-O')
                end
                hold off
                subplot(2, self.nch, k+self.nch)
                plot(self.Vm, cell2mat(self.ch_tau(k,:)))
                hold on;
                plot([-70, -70], [0 max(cell2mat(self.ch_tau(k,:)'))+.2],'r--')
                grid on; axis tight
                try
                legend('C-O','O-C')
                catch
                legend('C-O')
                end
                xlabel(self.channame{k})
                hold off
            end
            subplot(2, self.nch, 1)
            ylabel('eq state ratio')
            subplot(2, self.nch, self.nch+1)
            ylabel('act. t_const (ms)')
        end
    end
end