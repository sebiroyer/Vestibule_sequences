%function [trialL_N,serial_N,SxG_N,SxV_N] = MarkovModel(max_n_seg,n_mice,n_trials,start_positions,ini_RScSaS,P_transition)
%
%max_n_seg : the maximum number of segments in a trial (stop when the length of vestibule sequences exceed this number)
%Example: MarkovModel(10000,100,10,randi(23,10,1),[75 0 0 25],[50 20 5 25;10 80 0 10;10 10 50 30;10 30 10 50]);

function [trialL_N,serial_N,SxG_N,SxV_N] = MarkovModel(max_n_seg,n_mice,n_trials,start_positions,ini_RScSaS,P_transition)

%transition probabilities
R_RScSaS = P_transition(1,:);   %Random --> Random, Serial cw, Serial ccw, spatial
Sc_RScSaS = P_transition(2,:);  %Serial cw --> Random, Serial cw, Serial ccw, spatial
Sa_RScSaS = P_transition(3,:);  %Serial ccw --> Random, Serial cw, Serial ccw, spatial
S_RScSaS = P_transition(4,:);   %spatial --> Random, Serial cw, Serial ccw, spatial

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

        while e ~= goaldoorID && seg_ndx <= max_n_seg

            door_option = 0:23;
            toss_segment = randi([0 100],1);

            %Select strategy
            if seg_ndx == 1

                toss_segment = randi([0 100],1);
                if toss_segment <= ini_RScSaS(1)
                    strategy = 1;
                elseif toss_segment > ini_RScSaS(1) && toss_segment <= ini_RScSaS(1) + ini_RScSaS(2)
                    strategy = 2;
                elseif toss_segment > ini_RScSaS(1) + ini_RScSaS(2) && toss_segment <= ini_RScSaS(1) + ini_RScSaS(2) + ini_RScSaS(3)
                    strategy = 3;
                else
                    strategy = 4;
                end

            else

                if strategy == 1; P_RScSaS = R_RScSaS;
                elseif strategy == 2; P_RScSaS = Sc_RScSaS;
                elseif strategy == 3; P_RScSaS = Sa_RScSaS;
                else P_RScSaS = S_RScSaS;
                end

                if toss_segment <= P_RScSaS(1)
                    strategy = 1;
                elseif toss_segment > P_RScSaS(1) && toss_segment <= P_RScSaS(1) + P_RScSaS(2)
                    strategy = 2;
                elseif toss_segment > P_RScSaS(1) + P_RScSaS(2) && toss_segment <= P_RScSaS(1) + P_RScSaS(2) + P_RScSaS(3)
                    strategy = 3;
                else
                    strategy = 4;
                end

            end
            
            %Draw vestibule 
            if strategy == 1

                %1) random
                e = door_option(randi(24,1));

            elseif strategy == 2

                %2) serial clockwise
                mean_step1 = 1.2;
                std_step1 = 1.2;

                step = 0;
                while step == 0
                    step = round(mean_step1 + std_step1*randn);
                    if step<0; step = 0; end
                end
                e = b+step;
                if e>23
                    e = e-24;
                elseif e<0
                    e = 24+e;
                end

            elseif strategy == 3

                %3) serial counterclockwise
                mean_step2 = -2;
                std_step2 = 1.5;

                step = 0;
                while step == 0
                    step = round(mean_step2 + std_step2*randn);
                    if step>0; step = 0; end
                end
                e = b+step;
                if e>23
                    e = e-24;
                elseif e<0
                    e = 24+e;
                end

            else

                %4) spatial
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            MT_sbe = [MT_sbe;[mice_ii trial_ii seg_ndx b e]];
            seg_ndx = seg_ndx+1;
            b = e;

        end
    end
end


%compute distribution
n_seg = 10;

%gaps per segment_index
gaps = MT_sbe(:,5) - MT_sbe(:,4);
gaps(gaps>12) = -(24-gaps(gaps>12));
gaps(gaps<-12) = 24 + gaps(gaps<-12);

SxG_N = zeros(n_seg,25);
for ii = 1:n_seg
    sub_gaps = gaps(MT_sbe(:,3)==ii);
    n_seg_gap_tmp = histcounts(sub_gaps,-12.5:1:12.5);
    SxG_N(ii,:) = n_seg_gap_tmp/sum(n_seg_gap_tmp)*100;
end

%door visits per segment_index
SxV_N = zeros(n_seg,24);
for ii = 1:n_seg
    n_door_visits_tmp = histcounts(MT_sbe(MT_sbe(:,3)==ii,5),-0.5:1:23.5);
    n_door_visits_tmp = [n_door_visits_tmp(14:end) n_door_visits_tmp(1:13)];
    SxV_N(ii,:) = n_door_visits_tmp/sum(n_door_visits_tmp)*100;
end

%bout length and trial length
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

        %bout length
        A = cumsum([0;serial_move;0]);
        B = A(diff([0;serial_move;0])==-1);
        B = B - [0;B(1:end-1)];
        n = histcounts(B,0.5:1:20.5);
        TBl = [TBl;n];

        %n_seg_pertrial
        TNseg = [TNseg;length(ndx)];

    end

    MBl = [MBl;sum(TBl,1)];
    n = histcounts(TNseg,0.5:1:50.5);
    MNseg = [MNseg;n];
end
serial_N = mean(MBl,1);
serial_N = serial_N/sum(serial_N)*100;
trialL_N = mean(MNseg,1)/n_trials*100;


if nargout == 0
    figure;
    subplot(1,4,1);plot(1:20,serial_N,'b'); xlim([0 10]);xlabel('serial bout length');ylabel('% of serial bouts');
    subplot(1,4,2);plot(1:50,trialL_N,'b'); xlim([0 50]);xlabel('trial length (# of seg)');ylabel('% of trials');
    subplot(1,4,3);imagesc(-12:1:12,1:10,SxG_N); xlabel('seg length');ylabel('seg #');
    subplot(1,4,4);imagesc(-11:1:12,1:10,SxV_N); xlabel('door ID');ylabel('seg #');
end
