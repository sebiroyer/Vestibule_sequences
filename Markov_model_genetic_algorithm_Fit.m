%function Markov_model_genetic_algorithm_Fit(MTSBELDG,MTBl)
%
%

function Markov_model_genetic_algorithm_Fit(MTSBELDG,MTBl)

tic

%n_mice, n_trials, Day_range
Mice_range = 1:19;
n_mice = length(Mice_range);
Day_range = 6:15;       %after mice learned
n_trials = 10;
n_repetitions = 10;     %number of time to repeat the genetic algorithm simulation

%n_individuals thr_dist
n_individuals = 500;
n_generation = 500;
max_step = 10;          %maximum increment of probability values for mutations
thr_dist = 1;           %when generating new individuals, the difference between individuals should be > thr_dist for at least one of the probabilities
all_Pair_of_individuals = nchoosek(1:n_individuals/5,2);

%relative weigth of data points in the mean square error measure
mse_weigth = ones(1,50+20+10*25+10*24);       %trial length + serial bout + gap hist + visit hist
mse_weigth(1:10)=100;                           %first 10 bins of trial length
mse_weigth(50+(1:3))=100;                       %first 3 bins of serial bouts
mse_weigth(70+13+(0:25:225))=0;                 %peak of gap hist
mse_weigth(70+11+(0:25:225))=10;                %2 bins before 0 of gap hist
mse_weigth(70+14+(0:25:225))=10;                %1 bin after 0 of gap hist
mse_weigth(70+250+12+(0:24:216))=10;            %peaks of visit hist
mse_weigth(70+250+12)=500;                      %first peak of visit hist


%experimental data
[MTSBELDG,~,~,~,~,~,~,trialL_N,serial_N,SxG_N,SxV_N] = Distributions_segments(MTSBELDG,MTBl);
exp_conc = [trialL_N(:);serial_N(:);reshape(SxG_N',size(SxG_N(:)));reshape(SxV_N',size(SxV_N(:)))]';


%start positions
trial_range = zeros(1,length(Day_range)*n_trials);
for ii = 1:length(Day_range)
    trial_range((ii-1)*n_trials + (1:n_trials)) = Day_range(ii)*n_trials + (1:n_trials);
end
start_positions = MTSBELDG(MTSBELDG(:,1)==Mice_range(1) & ismember(MTSBELDG(:,2),trial_range) & MTSBELDG(:,3)==1,4);


%Genetic algorithm
for sim_nb = 1:n_repetitions

    %Generate probability matrices for the initial population of individuals
    P_transition = rand([5,4,n_individuals])-0.5;
    P_transition(P_transition<0) = 0.000001;
    P_transition = P_transition./repmat(sum(P_transition,2),1,4,1)*100;

    ndx = 10;
    while ~isempty(ndx)
        P_transition2D = reshape(P_transition,20,n_individuals);
        P_similarity = 1000*ones(n_individuals,n_individuals);
        for ii = 1:n_individuals-1
            for jj = ii+1:n_individuals
                P_similarity(ii,jj) = max(abs(P_transition2D(:,ii)-P_transition2D(:,jj)));
            end
        end
        P_similarity = P_similarity<thr_dist;
        ndx = find(sum(P_similarity,2)>0);

        if ~isempty(ndx)
            tmp = rand([5,4,length(ndx)]);
            tmp = tmp./repmat(sum(tmp,2),1,4,1)*100;
            P_transition(:,:,ndx) = tmp;
        end
    end


    %Make P_transition evolve: keep the best n_individuals/2 individuals and mutate/mixe them on each generation to generate an additional n_individuals/2 individuals
    all_P_transition2D = reshape(permute(P_transition,[2 1 3]),20,n_individuals)';
    all_best_mse = zeros(n_generation,n_individuals/2);
    all_best_P = zeros(n_generation,5,4,n_individuals/2);
    all_mse = zeros(n_individuals,1);
    for gg = 1:n_generation

        if gg == 1
            i_start = 1;
        else
            i_start = n_individuals/2+1;
        end

        for ii = i_start:n_individuals

            %run MarkovModel simulation and compute mean square error
            [trialL_Ns,serial_Ns,SxG_Ns,SxV_Ns] = MarkovModel(1000000,2,n_trials*length(Day_range),start_positions,P_transition(1,:,ii),P_transition(2:5,:,ii));
            conc = [trialL_Ns(:);serial_Ns(:);reshape(SxG_Ns',size(SxG_Ns(:)));reshape(SxV_Ns',size(SxV_Ns(:)))]';
            mse = mean(((exp_conc - conc).^2).*mse_weigth);

            %store
            all_mse(ii) = mse;

        end

        %keep the 50% best individuals and generate new individuals by mutating/mixing the 50% bests
        [all_mse,isort] = sort(all_mse);
        P_transition = P_transition(:,:,isort);
        Pair_of_individuals = all_Pair_of_individuals(randperm(size(all_Pair_of_individuals,1)),:);
        ndx = 1;
        for ii = 1:n_individuals/2
            D = 1;
            nloops = 0;
            while D~=0 && nloops<10
                nloops = nloops+1;

                if rem(gg,2)==0
                    %mutations: random increments (equivalent to gene mutation)
                    tmp = squeeze(P_transition(:,:,ii) + (rand([5,4])*(2*max_step)-max_step)); %randomly change all (works better)
                else
                    %mixing: combine parameters from pairs (equivalent to gene mixing)
                    ndx_para = unique(randi(5,4,1));
                    tmp = squeeze(P_transition(:,:,Pair_of_individuals(ndx,1)));
                    tmp(ndx_para,:) = squeeze(P_transition(ndx_para,:,Pair_of_individuals(ndx,2)));

                    ndx = ndx + 1;

                end

                %scale and store
                tmp(tmp<0)=0;
                tmp = tmp./repmat(sum(tmp,2),1,4)*100;
                tmp2D = reshape(tmp',1,20);
                D = sum(max(abs(all_P_transition2D - repmat(tmp2D,size(all_P_transition2D,1),1)),[],2)<thr_dist);
                if D==0 || nloops==10
                    P_transition(:,:,ii+n_individuals/2) = tmp;
                    all_P_transition2D = [all_P_transition2D;tmp2D];
                end
            end
        end

        all_best_mse(gg,:) = all_mse(1:n_individuals/2);
        all_best_P(gg,:,:,:) = P_transition(:,:,1:n_individuals/2);

    end

    if n_repetitions==1

        %generate simulation using the average of the 50% best P_transition from the last generation
        [trialL_Ns,serial_Ns,SxG_Ns,SxV_Ns] = MarkovModel(1000000,n_mice,n_trials*length(Day_range),start_positions,mean(P_transition(1,:,1:n_individuals/2),3),mean(P_transition(2:5,:,1:n_individuals/2),3));

        %plots experimental and simulation distributions
        figure;
        subplot(2,4,3);plot(1:20,serial_N,'b',1:20,serial_Ns,'r'); xlim([0 10]);xlabel('serial bout length');ylabel('% of serial bouts');
        subplot(2,4,4);plot(1:50,trialL_N,'b',1:50,trialL_Ns,'r'); xlim([0 50]);xlabel('trial length (# of seg)');ylabel('% of trials');
        subplot(2,4,5);imagesc(-12:1:12,1:10,SxG_N); xlabel('seg length');ylabel('seg #');title('exp')
        subplot(2,4,6);imagesc(-12:1:12,1:10,SxG_Ns); xlabel('seg length');ylabel('seg #');title('sim')
        subplot(2,4,7);imagesc(-11:1:12,1:10,SxV_N); xlabel('door ID');ylabel('seg #');title('exp')
        subplot(2,4,8);imagesc(-11:1:12,1:10,SxV_Ns); xlabel('door ID');ylabel('seg #');title('sim')

        figure;
        subplot(6,4,2:3);plot(1:n_generation,all_best_mse(:,1),'b',1:n_generation,mean(all_best_mse,2),'r');
        for ii = 1:5
            for jj = 1:4
                ndx = jj + ii*4;
                subplot(6,4,ndx);plot(1:n_generation,all_best_P(:,ii,jj,1),'b',1:n_generation,mean(all_best_P(:,ii,jj,:),4),'r');
                ylim([0 100]);
            end
        end

    end

    simulation_duration = toc;


    %save results
    save(['Markov_model_Fit_data' num2str(sim_nb) '.mat'],'all_best_P','all_best_mse','simulation_duration')

end
