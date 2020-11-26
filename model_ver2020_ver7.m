% model_ver2020_ver7.m:
% Created by Eugene M.Izhikevich.   2004
% Modified by Masanori Shiomono,    2012/2019
% Modified by Motoki Kajiwara,      2018
%
% http://www.izhikevich.org/publications/dastdp.htm
%
% This program reproduces the experiment in Fig.1 in
% Izhikevich E.M. (2007) Solving the Distal Reward Problem through Linkage
% of STDP and Dopamine Signaling. Cerebral Cortex, 10.1093/cercor/bhl152
%
% This model also involve inhibotory STDP.
% Wang, L., and Maffei, A. (2014).
% Inhibitory plasticity dictates the sign of plasticity at excitatory synapses.
% J. Neurosci. 34, 1083?1093. doi: 10.1523/JNEUROSCI.4711-13.2014
%
% n1 - the observing presynaptic neuron. syn is the synapse.
% Plot: top - spike raster. Bottom left - synaptic strength (blue), the
% eligibility trace (green), and the rewards (red x). Bottom right - the
% distribution of synaptic weights with the chosen synapse marked by red dot.

M=100;                 % number of synapses per neuron
D=4;                   % maximal conduction delay  (Originally 1ms, it's strange.)
% excitatory neurons   % inhibitory neurons      % total number
Ne=400;                Ni=100;                   N=Ne+Ni;
a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];
d=[   8*ones(Ne,1);    2*ones(Ni,1)];
smp =  15      % maximal positive synaptic strength
smn = -15      % maximal negative synaptic strength

post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]);

s=[5*ones(Ne,M);-4*ones(Ni,M)];     % synaptic weights
sd=zeros(N,M);                      % their derivatives

% parpool;%

for i=1:N
    if i<=Ne
        for j=1:D
            delays{i,j}=M/D*(j-1)+(1:M/D);
        end;
    else
        delays{i,1}=1:M;
    end;
    pre{i}=find(post==i & s>0);    % pre excitatory neurons
    aux{i}=N*(D-1-ceil(ceil(pre{i}/N)/(M/D)))+1+mod(pre{i}-1,N);
end;

STDP = zeros(N,1001+D);
v = -65*ones(N,1);                      % initial values
u = 0.2.*v;                             % initial values
firings=[-D 0];                         % spike timings

%% ---------------
% They are just for activity monitoring

T0              = 0.5*60*60*1;       % the starting time of data recording  (STDP ON   0.5 hour)
T_recstart = T0 + 30*60*1;           % the starting time of data recording  (STDP ON   0.5 hour)
T                = 1.5*60*60*1 + T0; % the duration of experiment           (STEP OFF  1.5 hour)


n1  = 1;             % a presynaptic excitatory neuron
n1h = 801;           % a presynaptic inhibitory neuron

syn  = 1;           % the synapse number to the exc. postsynaptic neuron
synh = 1;           % the synapse number to the inh. postsynaptic neuron

shist  = zeros(1000*T,2);
shisth = zeros(1000*T,2);

%% -----------------------

counter = ones(N,1);
totalfirings = zeros(N, (T-T_recstart)*1000);

% for long time sim.
% S = spalloc(m,n,nz)

disp('start!');
tic

for sec=1:T 
    
    if mod(sec, 50) == 0
        time_c = toc;
        
        if sec == 50
            disp(['Comp. time: ',num2str(T*(time_c/(60*sec))),'[min]']);
        end
        
        disp(['simulation time: ',num2str(sec/60),'[min]']);
        disp(['real spent time: ',num2str(time_c/60),'[min]']);
    end
    
    for t=1:1000                            % simulation of 1 sec
        % noise
        I=zeros(N,1);
        I(ceil(N*rand))=20; %
        
        fired = find(v>=30);                % indices of fired neurons
        v(fired)=-65;
        u(fired)=u(fired)+d(fired);
        STDP(fired,t+D)=0.1;
        
        for k=1:length(fired)
            sd(pre{fired(k)})=sd(pre{fired(k)})+STDP(N*t+aux{fired(k)});
        end
        
        %%
        if sec > T_recstart
            if sum(counter(fired) == length(totalfirings)) ~= 0
                totalfirings = [totalfirings zeros(N, (T-T_recstart)*1000)];
            end
            
            totalfirings(fired + counter(fired)*N) = t+1000*((sec-T_recstart)-1);
            counter(fired) = counter(fired) + 1;
        end
        
        firings=[firings;t*ones(length(fired),1),fired];
        k=size(firings,1);
        
        while firings(k,1)>t-D
            del=delays{firings(k,2),t-firings(k,1)+1};
            ind = post(firings(k,2),del);
            I(ind)=I(ind)+s(firings(k,2), del)';
            sd(firings(k,2),del)=sd(firings(k,2),del)-1.2*STDP(ind,t+D)';
            k=k-1;
        end;
        
        v=v+0.5*((0.04*v+5).*v+140-u+I);    % for numerical
        v=v+0.5*((0.04*v+5).*v+140-u+I);    % stability time
        u=u+a.*(0.2*v-u);                   % step is 0.5 ms
        STDP(:,t+D+1)=0.95*STDP(:,t+D);     % tau = 20 ms
        
        shist((sec-1)*1000+t,:)  = [mean(s(1:Ne,syn),1) , mean(sd(1:Ne,syn),1)];
        shisth((sec-1)*1000+t,:) = [mean(s(Ne+1:N,synh),1) , mean(sd(Ne+1:N,synh),1)];
    end
    
    %% ---- plot -------
    if mod(sec, 50) == 0
        figure(10+itr)
        subplot(2,1,1)
        plot(firings(:,1),firings(:,2),'.');
        axis([0 1000 0 N]);
        title(['Data ',num2str(data_index)],'fontsize', 15,'fontname','Arial');
        subplot(4,2,5);
        plot(0.001*(1:(sec-1)*1000+t),shist(1:(sec-1)*1000+t,1) );
        title(['Synaptic weight change (exc)'],'fontsize', 8,'fontname','Arial'); hold off;
        subplot(4,2,7);
        plot(0.001*(1:(sec-1)*1000+t),shisth(1:(sec-1)*1000+t,1));
        title(['Synaptic weight change (imh)'],'fontsize', 8,'fontname','Arial'); hold off;
        subplot(2,2,4);
        hist(s(find(s~=0)),2*smn:0.01:2*smp); % all synapses
        hold on; plot(s(n1,syn),0,'r.'); hold off;
        h = text(-3,5000,sprintf('firing rate: %f Hz', length(firings)/N));
        title(['All synaptic weight hist.'],'fontsize', 12,'fontname','Arial');
        set(gca,'XLim',[-12.5 12.5],'fontsize',12,'fontname','Arial');
        drawnow;
        delete(h);
    end
    
    %% ---- end plot ------
    
    STDP(:,1:D+1)=STDP(:,1001:1001+D);
    
    ind = find(firings(:,1) > 1001-D);
    firings=[-D 0;firings(ind,1)-1000,firings(ind,2)];
    
   if (mod(sec,10)==0)
        if sec <= T0
            s(1:Ne,:)   = max(0, min(smp,s(1:Ne,:)   + 0.002*sd(1:Ne,:)));     % exc. STDP, polycode shinya
            s(Ne+1:N,:) = min(0, max(smn,s(Ne+1:N,:) + 0.002*sd(Ne+1:N,:)));   % inh. STDP, polycode shinya
        end
        sd=0.9*sd; % this part is also very important and should locate out of "if sec <= T0" !!
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spikes = cell(N+2, 1);
for i=1:N
    spikes{i} =  full(totalfirings(i, find(totalfirings(i,:))));
end
spikes{end-1} = 1;
spikes{end} = [N, sec*1000];

save ./spikes spikes
