%function plot_Markov_model_Fit(MTSBELDG,MTBl)


function plot_Markov_model_Fit(MTSBELDG,MTBl)

n_sim = 10;

%Store mean values of all simulation repetitions 
all_mean_P = zeros(5,4,n_sim);
all_mean_mse = [];
for ss = 1:n_sim

    load(['Markov_model_Fit_data' num2str(ss) '.mat'])

    n_seeds = size(all_best_P,4);
    n_generation = size(all_best_P,1);

    for ii = 1:5
        for jj = 1:4
            all_mean_P(ii,jj,ss) = mean(all_best_P(n_generation,ii,jj,1:n_seeds/2),4);
        end
    end

    all_mean_mse = [all_mean_mse mean(all_best_mse,2)];

end
[~,i_Lowest_mean_mse] = min(all_mean_mse(n_generation,:));


%Plot probabilities for initial strategy and strategy transitions, for all repetitions (blue) and for best fit (red)
figure;
xx = repmat(1:4,1,n_sim);
for ii = 1:5
    subplot(5,1,ii);plot(xx,reshape(all_mean_P(ii,:,:),1,n_sim*4),'bo',(1:4)+0.1,all_mean_P(ii,:,i_Lowest_mean_mse),'ro')
    ylim([0 100])
    xlim([0 5])
    if ii == 1; ylabel('probability'); end
    if ii == 2; ylabel('R'); end
    if ii == 3; ylabel('Scw'); end
    if ii == 4; ylabel('Sccw'); end
    if ii == 5; ylabel('S'); xlabel('to R Scw Sccw S'); end
end

%Plot mean square error across generation
figure;
xx = 1:n_generation;
plot(xx,all_mean_mse,'k-')
ylabel('m.s.e'); xlabel('generation #');


%Generate Markov model distribution for the best fit average probabilities + plot experimental and simulated distributions
%Compute experimental distribution 
Mice_range = 1:19;
n_mice = length(Mice_range);
Day_range = 6:15;
n_trials = 10;
[~,~,~,~,~,~,~,trialL_N,serial_N,SxG_N,SxV_N] = Distributions_segments(MTSBELDG,MTBl);

%Compute start_positions
trial_range = zeros(1,length(Day_range)*n_trials);
for ii = 1:length(Day_range)
    trial_range((ii-1)*n_trials + (1:n_trials)) = Day_range(ii)*n_trials + (1:n_trials);
end
start_positions = MTSBELDG(MTSBELDG(:,1)==Mice_range(1) & ismember(MTSBELDG(:,2),trial_range) & MTSBELDG(:,3)==1,4);

%Run Markov model using the average Probabilities from the best fit
P = squeeze(all_mean_P(:,:,i_Lowest_mean_mse));
[trialL_Ns,serial_Ns,SxG_Ns,SxV_Ns] = MarkovModel(1000000,n_mice,n_trials*length(Day_range),start_positions,P(1,:),P(2:5,:));

%Plots experimental and simulated distributions
figure;
subplot(3,2,1);plot(1:20,serial_N,'b',1:20,serial_Ns,'r'); xlim([0 10]);xlabel('serial bout length');ylabel('% of serial bouts');
subplot(3,2,2);plot(1:50,trialL_N,'b',1:50,trialL_Ns,'r'); xlim([0 50]);xlabel('trial length (# of seg)');ylabel('% of trials');
subplot(3,2,3);imagesc(-12:1:12,1:10,SxG_N,[0 30]); xlabel('seg length');ylabel('seg #');title('exp')
subplot(3,2,4);imagesc(-12:1:12,1:10,SxG_Ns,[0 30]); xlabel('seg length');ylabel('seg #');title('sim')
subplot(3,2,5);imagesc(-11:1:12,1:10,SxV_N,[0 15]); xlabel('door ID');ylabel('seg #');title('exp')
subplot(3,2,6);imagesc(-11:1:12,1:10,SxV_Ns,[0 15]); xlabel('door ID');ylabel('seg #');title('sim')

