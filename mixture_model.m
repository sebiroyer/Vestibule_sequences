%function [n_seg_gap,n_door_visits,n_bouts,n_segDist] = mixture_model(n_mice,n_trials,start_positions,Percent_random,Percent_serial,n_conseq)
%
%Code for the mixture model
%Example: mixture_model(20,10,randi(23,10,1),30,30,6);


function [n_seg_gap,n_door_visits,n_bouts,n_segDist] = mixture_model(n_mice,n_trials,start_positions,Percent_random,Percent_serial,n_conseq)

%Goal ID, door ID, goal_door_distance
goaldoorID = 0;
doorID = 0:23;
goal_D = abs(doorID-goaldoorID);
goal_D12 = goal_D;
goal_D12(goal_D12<12) = -goal_D12(goal_D12<12);
goal_D12(goal_D12>12) = 24-goal_D12(goal_D12>12);


%simulation
MT_sbe = [];    %MT_sbe: [Mice_index Trial_index Segment_index Segment_start_position Segment_end_position]
for mice_ii = 1:n_mice
    for trial_ii = 1:n_trials

        b = start_positions(trial_ii);
        seg_ndx = 1;
        e = 100;

        while e ~= goaldoorID

            if mod(seg_ndx + n_conseq - 1,n_conseq) == 0
                toss_segment = randi([0 100],1);
            end
            door_option = 0:23;

            if toss_segment <= Percent_random

                %1) random strategy
                e = door_option(randi(24,1));

            elseif toss_segment <= Percent_random + Percent_serial

                %2) serial strategy
                mean_step1 = 1.2;
                std_step1 = 1.2;
                mean_step2 = -2;
                std_step2 = 1.5;
                percent_step_1 = 80;
                toss_steps = randi([0 100],1);

                step = 0;
                while step == 0
                    if toss_steps <= percent_step_1
                        step = round(mean_step1 + std_step1*randn);
                        if step<0; step = 0; end
                    else
                        step = round(mean_step2 + std_step2*randn);
                        if step>0; step = 0; end
                    end
                end
                e = b+step;
                if e>23
                    e = e-24;
                elseif e<0
                    e = 24+e;
                end

            else

                %3) spatial strategy
                reloop = 0;
                while reloop == 0
                    [door_option,i_d] = setdiff(0:23,b);    %no coming back to initial
                    goal_D12option = goal_D12(i_d);
                    mean_D = 0;
                    tau = 2;

                    goal_DD = round(mean_D + randi_expn(tau,1));

                    if goal_DD < -11
                        goal_DD = 24 + goal_DD;if goal_DD < -11;goal_DD = 24 + goal_DD;end
                    elseif goal_DD > 12
                        goal_DD = -(24 - goal_DD);if goal_DD > 12;goal_DD = -(24 - goal_DD);end
                    end
                    if sum(goal_D12option == goal_DD)>0
                        e = door_option(goal_D12option == goal_DD);
                        reloop = 1;
                    end
                end
            end

            MT_sbe = [MT_sbe;[mice_ii trial_ii seg_ndx b e]];
            seg_ndx = seg_ndx+1;
            b = e;

        end
    end
end


%compute distribution
gaps = MT_sbe(:,5) - MT_sbe(:,4);
gaps(gaps>12) = -(24-gaps(gaps>12));
gaps(gaps<-12) = 24 + gaps(gaps<-12);

n_seg_gap = histcounts(gaps,-12.5:1:12.5);
n_seg_gap = n_seg_gap/sum(n_seg_gap)*100;

n_door_visits = histcounts(MT_sbe(:,5),-0.5:1:23.5);
n_door_visits = [n_door_visits(14:end) n_door_visits(1:13)];
n_door_visits = n_door_visits/sum(n_door_visits)*100;

MBl = [];
MNseg = [];
for mice_ii = 1:n_mice
    TBl = [];
    TNseg = [];
    for trial_ii = 1:n_trials

        ndx = find(MT_sbe(:,1) == mice_ii & MT_sbe(:,2) == trial_ii);
        dd = MT_sbe(ndx,5)-MT_sbe(ndx,4);
        gg = abs(dd);
        gg(gg>=23)=1;
        serial_move = (gg == 1);  %serial move --> 1 door interspacing

        %serial bout length
        A = cumsum([0;serial_move;0]);
        B = A(diff([0;serial_move;0])==-1);
        B = B - [0;B(1:end-1)];
        n = histcounts(B,0.5:1:20.5);
        TBl = [TBl;n];

        %number of segments per trial (trial length)
        TNseg = [TNseg;length(ndx)];

    end

    MBl = [MBl;sum(TBl,1)];
    n = histcounts(TNseg,0.5:1:50.5);
    MNseg = [MNseg;n];
end
n_bouts = mean(MBl,1);
n_bouts = n_bouts/sum(n_bouts)*100;
n_segDist = mean(MNseg,1)/n_trials*100;


%Plot distributions if no output_argument specified 
if nargout == 0
    figure;
    subplot(1,4,1);plot(-12:1:12,n_seg_gap,'b'); xlim([-13 13]); ylim([0 max(n_seg_gap)*1.2]);xlabel('seg length');ylabel('% of segments');
    subplot(1,4,2);plot(-11:1:12,n_door_visits,'b'); xlim([-12 13]); ylim([0 max(n_door_visits)*1.2]);xlabel('door ID');ylabel('% of visits');
    subplot(1,4,3);plot(1:20,n_bouts,'b'); xlim([0 10]);xlabel('serial bout length');ylabel('% of serial bouts');
    subplot(1,4,4);plot(1:50,n_segDist,'b'); xlim([0 50]);xlabel('trial length (# of seg)');ylabel('% of trials');
end
