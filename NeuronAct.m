classdef NeuronAct < handle
    properties (SetAccess = protected)
        Vm;
        Ia;
        a;
        Im;
        fchan;
        gm;
        g;
        A;
        dt;
        Level;
        ds;
        Ca1;
        Ca2;
        ICa;
        nch;
        nf;
        ECa;
        conn;
        iRam;
        Cm;
        gs;
        Is;
        ram;
        ReshapeFun;
        Ic;
        Ra;
        Z;
        ncomp;
        spk;
        % channelstype;
    end
    methods
        function self=NeuronAct(channelstype,NeuronParams,ncell,dt)
            % initiation 
            % save for every 10^4 sample time.
            % I
            % self.channelstype=channelstype;
            self.nch=length(channelstype);
            self.ncomp=NeuronParams.numCompartments;
            nall=ncell*self.ncomp;
            if exist('restingV.mat','file')
                load('restingV.mat')
                self.Vm=V;
            else
                self.Vm=NeuronParams.E_leak*ones(nall,1)+rand(nall,1);
            end
            
            self.ReshapeFun=@(x)(reshape(x,ncell,self.ncomp));
            self.a=cell(self.nch,2);
            self.Ia=zeros(nall,1);
            self.Im=zeros(nall,self.nch);
            self.A=NeuronParams.A(:);
            self.dt=dt;
            self.Cm=NeuronParams.Cm(:);
            self.Level=2*ones(1,NeuronParams.numCompartments);
            self.Level(NeuronParams.compartmentParentArr==1)=1;
            self.Level(1)=0;
            self.Level=reshape(repmat(self.Level,ncell,1),[],1);
            [self.conn, self.ram]=ConnParts(NeuronParams);% the function of computing mean of V of nearby channels.
            self.Ra= 2./NeuronParams.Ra(:);
            self.iRam= 2./NeuronParams.Ra(:)+1./NeuronParams.Rm(:);
            self.Ic=NeuronParams.E_leak./NeuronParams.Rm(:);
            channames={'sdNa';...% Somato-dendritic Na+ channel (Na_7)
                %                 'aNa';...% Axonic Na+ channel (NaA_2)
                'wCa';...% High thereshold Ca++ channel (Ca_W)
                'tCa';...% Low thereshold Ca++ channel (CaT_3)
                'aKp';...% A type K+ channel, proximal (K_A_11)
                'aKd';...% A type K+ channel, distal (K_A_18)
                'drK';...% DR type K+ channel (K_DR_2)
                %                 'drKa';...% DR type K+ channel, axonic (K_DRA_4)
                'cK';...% C type K+ channel (K_C_1) ca
                'mK';...% M type K+ channel (K_M_4)
                'AHPK';...% AHP type K+ channel (K_AHP_Wtn)ca
                'Ih';...%  h type channel (Ih_3)
                };
            self.gm=[];%zeros(nall,self.nch);
            self.g=zeros(nall,self.nch);%cell(self.nch,1);
            self.ds=NeuronParams.ds(:);
            self.gs=0;
            self.Is=0;
            self.Z=reshape(repmat(mean(NeuronParams.compartmentZPositionMat,2)',ncell,1),[],1);
            
            for k=1:self.nch
                channame=channames{channelstype(k)};
                self.fchan{k}=str2func(channame);
                self=self.fchan{k}(self,k);
            end
            w= 10^-11 *((NeuronParams.compartmentDiameterArr>2)+((NeuronParams.compartmentDiameterArr<=2).*NeuronParams.compartmentDiameterArr)); % cm, thickness of the diffusion shell
            self.nf=1/2/96.48533./bsxfun(@times,NeuronParams.A,w');%  Ca influx affecting the i-th pool
            self.nf=self.nf(:);
            self=updateCalsium(self,0);
            self.spk=false(ncell,1);
        end
        %% update channels.
        function self=updateNeuron(self,Gglu,GgabaA,GgabaB)
            %  self=updateNeuron(self,Gglu,GgabaA,GgabaB)
            self.Is=-reshape(self.ram.*(self.ReshapeFun(self.Vm.*self.Ra)*self.conn),[],1)-self.Ic;
            self.gs=self.iRam;
            for k=1:self.nch
                
                self=self.fchan{k}(self,k);
                if sum(isnan(self.gs)|isinf(self.gs))>0
                    fprintf('wrong at %d',k)
                    break
                end
            end
            self=updateCalsium(self,1);
            
            % --- input current --- %
            self.Is=self.Is+GgabaA*10.71+GgabaB*3;%75/7 90/30
            self.gs=self.gs+GgabaA/7+GgabaB/30+Gglu/2;
            %             IgabaA=Inputs.GgabaA.*(self.Vm+75);
            %             IgabaB=Inputs.GgabaB.*(self.Vm+90);
            %             Iglu=Inputs.Gglu.*self.Vm;
            % -- cable conductance --%
            % [ga,gI]=self.conn(self.Vm,self.Ra);
            % self.Is=self.Is+self.ram.*(self.Vm*self.conn);
            % self.gs=self.gs-self.iRa; % constant!
            self.Vm=self.Vm-self.dt*(self.Is+self.gs.*self.Vm)./self.Cm;
            tm=self.Vm(1:self.ncomp:end);
            %              self.spk=(tm>-40);
            %              tm(self.spk)=-70;
            self.Vm(1:self.ncomp:end)=tm;
            
            
            %             self.Vm=(self.Vm+(self.Is./self.gs)).*exp(-self.dt*self.gs./self.Cm)-(self.Is./selfgs);
        end
        %% Intracellular Ca++ dynamics
        function self=updateCalsium(self,k)
            if k==0
                self.Ca1=0.05;%0;%
                self.Ca2=0.05;%0;%
                self.ICa=0;
                self.ECa=-13.275*log((max(self.Ca1,10^-9))/1200);%  +.05[C
            else
                % FA=96.48533;% mC/ uM
                % z=2;
                % tau=[.9; 10^3] ;% ms,  time constant
                % w= 10^-3 *((NP.diam>2)+((NP.diam<2).*NP.diam)); % cm, thickness of the diffusion shell
                ab=(.05/.9-0.7*self.nf.*self.ICa )*.9;%/1000mean()
                self.Ca1=(self.Ca1-ab).*exp(-self.dt/.9)+ab;
                % ab=(.05/.9-self.nf.*self.ICa )*.9;%/1000mean()
                % self.Cai=(self.Cai-ab).*exp(-self.dt/.9)+ab;
                ab=(.05/10^3 -0.024*self.nf.*self.ICa)*10^3;%mean()
                self.Ca2=(self.Ca2-ab).*exp(-self.dt/10^3)+ab;
                % Cabase=.05;% uM
            end
            self.ECa=-13.275*log((max(self.Ca1,10^-9)+.05)/1200);% +.05 [Ca]o w
            self.ICa=self.ICa*0;
            % where taui is the removal time constant in the i-th calcium pool, fi the
            % fractio
        end
        %% ---------------Somato-dendritic Na+ channel (Na_7)-------
        function self=sdNa(self,k)
            if size(self.gm,2)<k
                self.gm(:,k)=self.A.*500/100^2;
                self.a{k,1}=0;
                self.a{k,2}=0;
                %             gNa=500/100^2;% S/cm2
            else
                minf = 1./(1+exp((25+self.Vm)/-5));
                taum = .020+.22./(exp((self.Vm+55)/13)+exp((self.Vm+55)/-80));% ms
                % dm=dt./taum.*(minf-m);
                self.a{k,1}=(self.a{k,1}-minf).*exp(-self.dt./taum)+minf;
                
                hinf = 1./(1+exp((50+self.Vm)/4));
                tauh = .6+8*exp((self.Vm+45)/13)./(1+exp((self.Vm+45)/5.5));
                %  dh=dt./tauh.*(hinf-h);
                self.a{k,2}=(self.a{k,2}-hinf).*exp(-self.dt./tauh)+hinf;
                self.g(:,k)=self.gm(:,k).*self.a{k,1}.^3.*self.a{k,2};
                self.Im(:,k)= self.g(:,k).*(self.Vm-50);% record if you like. Recording
                
                self.gs=self.gs+self.g(:,k);
                self.Is=self.Is-50*self.g(:,k);
            end
            % Density: 500 S/m2 in soma and 500 S/m2 all over the apical and basal tree
        end
        %% ---------------------Axonic Na+ channel (NaA_2)-------
        function self=aNa(self,k)
            if size(self.gm,2)<k
                self.gm(:,k)=self.A.*((self.Level==1)*666+(self.Level>1)*1000)/100^2;% S/cm2
                self.a{k,1}=0;
                self.a{k,2}=0;
            else
                minf = 1./(1+exp((32+self.Vm)/-4.2));
                taum = .020+.220./(exp((self.Vm+62)/13)+exp((self.Vm+62)/-80));
                % dm=dt./taum.*(minf-m);
                self.a{k,1}=(self.a{k,1}-minf).*exp(-self.dt./taum)+minf;
                
                
                hinf = 1./(1+exp((50+self.Vm)/4));
                tauh = .0200+8*exp((self.Vm+45)/13)./(1+exp((self.Vm+45)/5.5));
                %  dh=dt./tauh.*(hinf-h);
                self.a{k,2}=(self.a{k,2}-hinf).*exp(-self.dt/tauh)+hinf;
                self.g(:,k)=self.gm(:,k).*self.a{k,1}.^3.*self.a{k,2};
                self.Im(:,k)= self.g(:,k).*(self.Vm-50);% record if you like. Recording
                
                self.gs=self.gs+self.g(:,k);
                self.Is=self.Is-50*self.g(:,k);
            end
            % Density: 666 S/m2 in the nearest compartment to soma and 1000 S/m2 in the unmielinized compartments thereof
        end
        
        %% ----------------High thereshold Ca++ channel (Ca_W)--------------
        function self=wCa(self,k)
            if size(self.gm,2)<k
                self.gm(:,k)=self.A.*((self.ds<=100)|self.Z>0).*30/100^2;% S/cm2
                self.a{k,1}=0;
                self.a{k,2}=0;
            else
                self.Vm(abs(self.Vm-26)<10^-4 |abs(self.Vm-12)<10^-4 )=self.Vm(abs(self.Vm-26)<10^-4 |abs(self.Vm-12)<10^-4 )+10^-3;
                alpha = .160*(self.Vm+26)./(-exp((self.Vm+26)/-4.5)+1);
                beta = .040*(self.Vm+12)./(exp((self.Vm+12)/10)-1);
                % ds=alphas.*(1-s)+betas.*s;
                minf=alpha./(alpha-beta);
                self.a{k,1}=(self.a{k,1}-minf).*exp(-self.dt.*(alpha-beta))+minf;
                
                
                alpha = 2./exp((self.Vm+94)/10);
                beta = 8./(exp((self.Vm-68)/-27)+1);
                % dr=alpha.*(1-r)+beta.*r;
                hinf=alpha./(alpha-beta);
                self.a{k,2}=(self.a{k,2}-hinf).*exp(-self.dt.*(alpha-beta))+hinf;
                self.g(:,k)=self.gm(:,k).*self.a{k,1}.^2.*self.a{k,2};
                ii=self.g(:,k).*(self.Vm-self.ECa);% r
                self.Im(:,k)= ii;%self.g(:,k).*(self.Vm-self.ECa);% record if you like. Recording
                
                self.ICa=self.ICa+ii;
                self.gs=self.gs+self.g(:,k);
                self.Is=self.Is-self.ECa.*self.g(:,k);
            end
            % Density: 666 S/m2 in the nearest compartment to soma and 1000 S/m2 in the unmielinized compartments thereof
        end
        %% --------------------------Low thereshold Ca++ channel (CaT_3)----
        function self=tCa(self,k)
            if size(self.gm,2)<k
                self.gm(:,k)=self.A.*min((self.Z<0).*(self.ds>100).*(30 + 25*(self.ds/100-1)),100)/100^2;%S/cm2
                self.a{k,1}=0;
                self.a{k,2}=0;
            else
                sinf = 1./(1+exp((54+self.Vm)/-5));
                taus = .204+.333./(exp((self.Vm+15.8)/18.2)+exp((self.Vm+131)/-16.7));
                self.a{k,1}=(self.a{k,1}-sinf).*exp(-self.dt./taus)+sinf;
                
                rinf = 1./(1+exp((80+self.Vm)/4));
                taur=rinf;
                taur(self.Vm<= -81)=.333*exp((self.Vm(self.Vm<= -81)+466)/66.6);
                taur(self.Vm>-81)=9.32 + .333*exp((self.Vm(self.Vm> -81)+21)/-10.5);
                
                self.a{k,2}=(self.a{k,2}-rinf).*exp(-self.dt./taur)+rinf;
                self.g(:,k)=self.gm(:,k).*self.a{k,1}.^2.*self.a{k,2};
                ii=self.g(:,k).*(self.Vm-self.ECa);% record if you like. Recording
                self.Im(:,k)= ii;
                self.ICa=self.ICa+ii;
                
                self.gs=self.gs+self.g(:,k);
                self.Is=self.Is-self.ECa.*self.g(:,k);
            end
            % Density: 666 S/m2 in the nearest compartment to soma and 1000 S/m2 in the unmielinized compartments thereof
        end
        %% -------------------A type K+ channel, proximal (K_A_11)-------------
        function self=aKp(self,k)
            if size(self.gm,2)<k
                self.gm(:,k)=self.A.*(self.ds<100).*((self.Z>0).*(20+2*self.ds)+(self.Z<0).*(20+1.7*self.ds))/100^2;%S/cm2
                self.a{k,1}=0;
                self.a{k,2}=0;
            else
                ninf = 1./(1+exp((12+self.Vm)/-8.5));
                taun = .080;
                self.a{k,1}=(self.a{k,1}-ninf).*exp(-self.dt./taun)+ninf;
                
                
                linf = 1./(1+exp((56+self.Vm)/7));
                taul = max(.26*(self.Vm+50),2);
                self.a{k,2}=(self.a{k,2}-linf).*exp(-self.dt./taul)+linf;
                self.g(:,k)=self.gm(:,k).*self.a{k,1}.*self.a{k,2};
                self.Im(:,k)= self.g(:,k).*(self.Vm+90);% record if you like. Recording
                self.gs=self.gs+self.g(:,k);
                self.Is=self.Is+90*self.g(:,k);
            end
            % Density: 666 S/m2 in the nearest compartment to soma and 1000 S/m2 in the unmielinized compartments thereof
        end
        
        %% ------------------ A type K+ channel, distal (K_A_18)---------------
        function self=aKd(self,k)
            if size(self.gm,2)<k
                self.gm(:,k)=self.A.*min((self.ds>100).*((self.Z>0).*(20+2*self.ds)+(self.Z<0).*(20+1.7*self.ds)),600)/100^2;%S/cm2
                self.a{k,1}=0;
                self.a{k,2}=0;
            else
                ninf = 10.^(self.Vm<= -65)./(1+exp((24+self.Vm)/-7));
                taun = .080;
                self.a{k,1}=(self.a{k,1}-ninf).*exp(-self.dt./taun)+ninf;
                
                linf = 1./(1+exp((56+self.Vm)/7));
                taul = max(0.26*(self.Vm+50),2);
                self.a{k,2}=(self.a{k,2}-linf).*exp(-self.dt./taul)+linf;
                self.g(:,k)=self.gm(:,k).*self.a{k,1}.*self.a{k,2};
                self.Im(:,k)= self.g(:,k).*(self.Vm+90);% record if you like. Recording
                self.gs=self.gs+self.g(:,k);
                self.Is=self.Is+90*self.g(:,k);
            end
            % Density: 666 S/m2 in the nearest compartment to soma and 1000 S/m2 in the unmielinized compartments thereof
        end
        
        %% ---------------DR type K+ channel (K_DR_2)-------------------
        function self=drK(self,k)
            if size(self.gm,2)<k
                self.gm(:,k)=self.A.*40/100^2;%S/cm2
                self.a{k,1}=0;
                self.a{k,2}=0;
            else
                ninf = 1./(1+exp((5-self.Vm)/11));
                taun = 1.2;
                self.a{k,1}=(self.a{k,1}-ninf).*exp(-self.dt./taun)+ninf;
                
                self.g(:,k)=self.gm(:,k).*self.a{k,1};
                self.Im(:,k)= self.g(:,k).*(self.Vm+90);% record if you like. Recording
                self.gs=self.gs+self.g(:,k);
                self.Is=self.Is+90*self.g(:,k);
            end
            % Density: 666 S/m2 in the nearest compartment to soma and 1000 S/m2 in the unmielinized compartments thereof
        end
        
        %% ----------------C type K+ channel (K_C_1)----------------
        function self=cK(self,k)
            if size(self.gm,2)<k
                self.gm(:,k)=self.A.*((self.Z>0 & self.ds<225).*(720+3.2*self.ds)+(self.Z==0)*720+(self.Z<0 & self.ds<100).*(720+7*self.ds))/100^2;% S/cm2
                self.a{k,1}=0;
                self.a{k,2}=0;
            else
                Vshift = 40.*log(max((self.Ca1+.05)*1000,10^-7))-105;
                alpha = -7.7*(self.Vm+Vshift+103)/1000./(exp((self.Vm+Vshift+103)/-12)-1);
                beta = 1.7./exp((self.Vm+Vshift+237)/30);
                % dc=alphac.*(1-c)+betac.*c;
                minf=alpha./(alpha-beta);
                self.a{k,1}=(self.a{k,1}-minf).*exp(-self.dt.*(alpha-beta))+minf;
                
                
                alpha = 1./exp((self.Vm+79)/10);
                beta =4./(exp((self.Vm-82)/-27)+1);
                % dd=alphad.*(1-d)+betad.*d;
                hinf=alpha./(alpha-beta);
                self.a{k,2}=(self.a{k,2}-hinf).*exp(-self.dt.*(alpha-beta))+hinf;
                self.g(:,k)=self.gm(:,k).*self.a{k,1}.^2.*self.a{k,2};
                self.Im(:,k)= self.g(:,k).*(self.Vm+90);% record if you like. Recording
                self.gs=self.gs+self.g(:,k);
                self.Is=self.Is+90*self.g(:,k);
            end
            % Density (ds = distance to soma in µm):
            % Soma: 720 S/m2
            % Apical tree:
            %     If ds < 225, Dens = 720 - 320(ds/100) S/m2
            %     If ds > 225, Dens = 0.0 S/m2
            % Basal tree:
            %     If ds < 100, Dens = 720 - 700(ds/100) S/m2
            %     If ds > 100, Dens = 0.0 S/m2
            
            
        end
        
        %% --------------- M type K+ channel (K_M_4)--------------------
        function self=mK(self,k)
            if size(self.gm,2)<k
                self.gm(:,k)=self.A.*((self.Z>0).*5+(self.Z==0)*10+(self.Z<0).*4.5)/100^2;%S/cm2
                self.a{k,1}=0;
                self.a{k,2}=0;
            else
                ninf = 1./(1+exp((40+self.Vm)/-8.5));
                taun =0.9+20*exp((self.Vm+38)/25)./(1+exp((self.Vm+38)/7));
                self.a{k,1}=(self.a{k,1}-ninf).*exp(-self.dt./taun)+ninf;
                self.g(:,k)=self.gm(:,k).*self.a{k,1};
                self.Im(:,k)= self.g(:,k).*(self.Vm+90);% record if you like. Recording
                self.gs=self.gs+self.g(:,k);
                self.Is=self.Is+90*self.g(:,k);
            end
            % Density: 666 S/m2 in the nearest compartment to soma and 1000 S/m2 in the unmielinized compartments thereof
        end
        %% ----------------AHP type K+ channel (K_AHP_Wtn)----------------
        function self=AHPK(self,k)
            if size(self.gm,2)<k
                self.gm(:,k)=self.A.*max((self.Z>0).*(6 -.05*self.ds)+(self.Z<0).*(6 -.08*self.ds),1)/100^2;% S/cm2
                self.a{k,1}=0;
                self.a{k,2}=0;
            else
                alpha = 4.8*10^6./exp(-5*log(max((self.Ca2+.05)*1000-35,10^-6)));
                beta = 12*10^6./exp(2*log(max((self.Ca2+.05)*1000+100,10^-6)));
                %dq=alphaq.*(1-q)+betaq.*q;
                minf=alpha./(alpha-beta);
                self.a{k,1}=(self.a{k,1}-minf).*exp(-self.dt.*(alpha-beta))+minf;
                
                self.g(:,k)=self.gm(:,k).*self.a{k,1};
                self.Im(:,k)= self.g(:,k).*(self.Vm+90);% record if you like. Recording
                self.gs=self.gs+self.g(:,k);
                self.Is=self.Is+90*self.g(:,k);
            end
            % Density (ds = distance to soma in µm):
            % Soma: 6 S/m2
            % Apical tree:
            %     Dens = 6 - 5(ds/100) S/m2
            % Basal tree:
            %     Dens = 6 - 8(ds/100) S/m2
            % If Dens < 1.0 S/m2, Dens = 1.0 S/m2
        end
        %% ---------------------h type channel (Ih_3)---------------------
        function self=Ih(self,k)
            if size(self.gm,2)<k
                self.gm(:,k)=(self.A.*(self.Z>0).*(4 + 116./(1.0 + exp((250 - self.Z)/100)))+...
                    self.A.*(self.Z<0).*(4 + 96./(1.0 + exp((250 - self.Z)/100)))+...
                    self.A.*(self.Z==0).*4)/100^2;%S/cm2
                
                self.a{k,1}=0;
                self.a{k,2}=0;
            else
                ninf = 1./(1+exp((90+self.Vm)/4));
                taun =0.6+40.*exp((self.Vm+65)/35)./(1+exp((self.Vm+65)/3.5));
                self.a{k,1}=(self.a{k,1}-ninf).*exp(-self.dt./taun)+ninf;
                self.g(:,k)=self.gm(:,k).*self.a{k,1};
                self.Im(:,k)= self.g(:,k).*(self.Vm+90);% record if you like. Recording
                self.gs=self.gs+self.g(:,k);
                self.Is=self.Is+90*self.g(:,k);
            end
            % Density: 666 S/m2 in the nearest compartment to soma and 1000 S/m2 in the unmielinized compartments thereof
        end
        
    end
end