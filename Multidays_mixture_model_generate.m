%function Multidays_mixture_model_generate(MTSBELDG)
%
%Code to generate simulated data for various values of P_random P_serial and P_spatial using the mixture model
%Note that P_spatial is note explicitely specified since P_spatial = 100 - P_random - P_serial
%Simply put the experimental data matrix MTSBELDG as input

function Multidays_mixture_model_generate(MTSBELDG)


%simulation parameters
n_mice = 20;
n_trials = 10;
Day_range = 1:19;
Percent_random = 0:2:100;
Percent_serial = 0:2:100;
Gx10xPRxPSxD = nan(25,10,length(Percent_random),length(Percent_serial),length(Day_range)); %segment length
NDVx10xPRxPSxD = nan(24,10,length(Percent_random),length(Percent_serial),length(Day_range)); %visited door ID
NSBx10xPRxPSxD = nan(20,10,length(Percent_random),length(Percent_serial),length(Day_range)); %serial bouts
Nsegx10xPRxPSxD = nan(50,10,length(Percent_random),length(Percent_serial),length(Day_range)); %trial length
n_segments_before_strategy_switch = 6;


%run mixture model simulation for each day, each values of Percent_random/Percent_serial, and 10 repetitions
for dd = 1:length(Day_range)

    %get start_positions for that day, from experimental data 
    trial_range = Day_range(dd)*n_trials + (1:n_trials);
    start_positions = MTSBELDG(MTSBELDG(:,1)==Mice_range(1) & ismember(MTSBELDG(:,2),trial_range) & MTSBELDG(:,3)==1,4);

    for ii = 1:length(Percent_random)
        for jj = 1:length(Percent_serial)
            if Percent_random(ii) + Percent_serial(jj) <=100
                for kk = 1:10
                    [n_seg_gap,n_door_visits,n_bouts,n_segDist] = mixture_model(n_mice,n_trials,start_positions,Percent_random(ii),Percent_serial(jj),n_segments_before_strategy_switch);
                    Gx10xPRxPSxD(:,kk,ii,jj,dd) = n_seg_gap;
                    NDVx10xPRxPSxD(:,kk,ii,jj,dd) = n_door_visits;
                    NSBx10xPRxPSxD(:,kk,ii,jj,dd) = n_bouts;
                    Nsegx10xPRxPSxD(:,kk,ii,jj,dd) = n_segDist;
                end
            end
        end
    end

end

%store simulated data matrices into a file 
save('simulated_data','Gx10xPRxPSxD','NDVx10xPRxPSxD','NSBx10xPRxPSxD','Nsegx10xPRxPSxD')