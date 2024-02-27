%function [MTSBELDG,GxD_N,VxD_N,NxD_serial,D_NsegDist,G_N,V_N,trialL_N,serial_N,SxG_N,SxV_N] = Distributions_segments(MTSBELDG,MTBl)
%
%

function [MTSBELDG,GxD_N,VxD_N,NxD_serial,D_NsegDist,G_N,V_N,trialL_N,serial_N,SxG_N,SxV_N] = Distributions_segments(MTSBELDG,MTBl)


n_trials = 10;                                                  %number of trials per day
n_days = length(unique(floor((MTBl(:,2)-1)/n_trials)));         %day1-->trial11:20, day2-->trial21:30
n_mice = length(unique(MTBl(:,1)));

day_range = 6:15;

segments_bins_edges = 0.5:1:10.5;   

%mattrices across segments
MxDxTxG_PL = nan(n_mice,n_days,10,25);                                          %mean segment path length, for each mouse, day, trial, gapsize
MxDxTxG_N = nan(n_mice,n_days,10,25);                                           %segment number, for each mouse, day, trial, gapsize
MxDxTxG_D = nan(n_mice,n_days,10,25);                                           %mean segment duration, for each mouse, day, trial, gapsize
MxDxSxG_N = nan(n_mice,n_days,length(segments_bins_edges)-1,25);                %segment number, for each mouse, day, segmentNb, gapsize
MxDxT_Nseg = nan(n_mice,n_days,10);                                             %mean segment number, for each mouse, day, trial
MxD_NsegDist = nan(n_mice,n_days,50);                                           %distribution of segment number, for each mouse,day
MxD_door_visits = zeros(n_mice,n_days,24);                                      %number of visits, for each mouse, day, door
MxDxS_door_visits = zeros(n_mice,n_days,length(segments_bins_edges)-1,24);      %number of visits, for each mouse, day, segmentNb, door
MxD_n_serial = zeros(n_mice,n_days,20);                                         %distribution of serial bout length, for each mouse, day

door_gaps = -12:1:12;
for ii = 1:n_days
    for mice_ii = 1:n_mice

        %number of segments by gap size, and path length by gap size
        for jj = 1:25
            for kk = 1:10
                ndx = find(MTSBELDG(:,8)==door_gaps(jj) & MTSBELDG(:,2)==kk+ii*10 & MTSBELDG(:,1)==mice_ii);
                MxDxTxG_N(mice_ii,ii,kk,jj) = length(ndx); %per trial
                if ~isempty(ndx)
                    MxDxTxG_PL(mice_ii,ii,kk,jj) = mean(MTSBELDG(ndx,6));
                    MxDxTxG_D(mice_ii,ii,kk,jj) = mean(MTSBELDG(ndx,7));
                end                
            end
            
            %number of segments by gap size, per segment bin
            for kk = 1:length(segments_bins_edges)-1
                ndx = find(MTSBELDG(:,3)>segments_bins_edges(kk) & MTSBELDG(:,3)<=segments_bins_edges(kk+1) & MTSBELDG(:,8)==door_gaps(jj) & MTSBELDG(:,2)>ii*10 & MTSBELDG(:,2)<=(ii+1)*10 & MTSBELDG(:,1)==mice_ii);
                MxDxSxG_N(mice_ii,ii,kk,jj) = length(ndx);
            end
        end

        %number of segment per trial
        for jj = 1:10
            ndx = find(MTSBELDG(:,2)==(ii*10+jj) & MTSBELDG(:,1)==mice_ii);
            MxDxT_Nseg(mice_ii,ii,jj) = length(ndx);
        end
        MxD_NsegDist(mice_ii,ii,:) = histcounts(MxDxT_Nseg(mice_ii,ii,:),0.5:1:50.5);
        
        %door visits across days
        ndx = find(MTSBELDG(:,1)==mice_ii & MTSBELDG(:,2)>ii*10 & MTSBELDG(:,2)<=(ii+1)*10);
        tmp = histcounts(MTSBELDG(ndx,5),-0.5:1:23.5);
        MxD_door_visits(mice_ii,ii,:) = [tmp(14:end) tmp(1:13)];

        %door visits across days, per segment bin 
        for kk = 1:length(segments_bins_edges)-1
            ndx = find(MTSBELDG(:,1)==mice_ii & MTSBELDG(:,2)>ii*10 & MTSBELDG(:,2)<=(ii+1)*10 & MTSBELDG(:,3)>segments_bins_edges(kk) & MTSBELDG(:,3)<=segments_bins_edges(kk+1));
            tmp = histcounts(MTSBELDG(ndx,5),-0.5:1:23.5);
            MxDxS_door_visits(mice_ii,ii,kk,:) = [tmp(14:end) tmp(1:13)];
        end

        %n_serial across days
        ndx = find(MTBl(:,1)==mice_ii & MTBl(:,2)>ii*10 & MTBl(:,2)<=(ii+1)*10);
        MxD_n_serial(mice_ii,ii,:) = sum(MTBl(ndx,3:end),1);
        
    end
end

%across days
GxD_PL = squeeze(mean(MxDxTxG_PL,[1 3],'omitnan'))';
GxD_D = squeeze(mean(MxDxTxG_D,[1 3],'omitnan'))';
GxD_N = squeeze(mean(MxDxTxG_N,[1 3],'omitnan'))';
GxD_N = GxD_N./repmat(sum(GxD_N,1),size(GxD_N,1),1)*100;
VxD_N = squeeze(mean(MxD_door_visits,1,'omitnan'))';
VxD_N = VxD_N./repmat(sum(VxD_N,1),size(VxD_N,1),1)*100;
MxD = squeeze(mean(MxDxT_Nseg,3,'omitnan'));
D_Nseg = mean(MxD,1);
D_NsegDist = squeeze(mean(MxD_NsegDist,1))'/10*100;
NxD_serial = squeeze(mean(MxD_n_serial(:,:,1:20),1))';
NxD_serial = NxD_serial./repmat(sum(NxD_serial,1),size(NxD_serial,1),1)*100;

%average over day_range
G_N = squeeze(mean(MxDxTxG_N(:,day_range,:,:),1:3,'omitnan'));
G_N = G_N/sum(G_N)*100;
V_N = squeeze(mean(MxD_door_visits(:,day_range,:),1:2,'omitnan'));
V_N = V_N/sum(V_N)*100;
trialL_N = squeeze(mean(MxD_NsegDist(:,day_range,:),1:2))/10*100;
serial_N = squeeze(mean(MxD_n_serial(:,day_range,1:20),1:2));
serial_N = serial_N/sum(serial_N)*100;

%across segments
SxG_N = squeeze(mean(MxDxSxG_N(:,day_range,:,:),[1 2],'omitnan'));
SxG_N = SxG_N./repmat(sum(SxG_N,2),1,size(SxG_N,2))*100;
SxV_N = squeeze(mean(MxDxS_door_visits(:,day_range,:,:),[1 2],'omitnan'));
SxV_N = SxV_N./repmat(sum(SxV_N,2),1,size(SxV_N,2))*100;


if nargout == 0
    D1 = 1;
    D2 = 6:15;
    figure;
    subplot(5,3,1);imagesc(GxD_PL);ylabel('seg length');xlabel('days');
    subplot(5,3,2);plot(-12:1:12,mean(GxD_PL(:,D1),2,'omitnan'),'b',-12:1:12,mean(GxD_PL(:,D2),2,'omitnan'),'r');ylabel('path length (cm)');xlabel('seg length');
    subplot(5,3,3);imagesc(GxD_D);ylabel('seg length');xlabel('days');
    subplot(5,3,4);imagesc(GxD_N);ylabel('seg length');xlabel('days');
    subplot(5,3,5);plot(-12:1:12,mean(GxD_N(:,D1),2),'b',-12:1:12,mean(GxD_N(:,D2),2),'r');ylabel('% of segments');xlabel('seg length');
    subplot(5,3,7);imagesc(VxD_N);ylabel('door#');xlabel('days');
    subplot(5,3,8);plot(-11:12,mean(VxD_N(:,D1),2)/sum(mean(VxD_N(:,D1),2))*100,'b',-11:12,mean(VxD_N(:,D2),2)/sum(mean(VxD_N(:,D2),2))*100,'r');ylabel('% of visits');xlabel('door#');
    subplot(5,3,10);imagesc(NxD_serial);xlabel('days');ylabel('serial bout length');
    subplot(5,3,11);plot(1:20,mean(NxD_serial(:,D1),2),'b',1:20,mean(NxD_serial(:,D2),2),'r');xlabel('serial bout length');ylabel('% of serial bouts');
    subplot(5,3,13);imagesc(D_NsegDist);xlabel('days');ylabel('trial length');
    subplot(5,3,14);plot(1:50,mean(D_NsegDist(:,D1),2),'b',1:50,mean(D_NsegDist(:,D2),2),'r');xlabel('trial length (seg#)');ylabel('# of trials');
    subplot(5,3,15);plot(1:15,MxD,'b.',1:15,mean(MxD,1),'ro-');ylabel('trial length (seg#)');xlabel('days');

    figure;
    subplot(3,2,1);imagesc(SxG_N);xlabel('seg length');ylabel('segment');
    subplot(3,2,3);plot(1:25,SxG_N(1,:),'b',1:25,SxG_N(2,:),'r',1:25,SxG_N(3,:),'g');ylabel('%');xlabel('seg length');
    subplot(3,2,5);plot(1:10,SxG_N(:,14),'b');ylim([0 50]);xlabel('% of segment size==1');ylabel('segment');
    subplot(3,2,2);imagesc(SxV_N);xlabel('door#');ylabel('segment');
    subplot(3,2,4);plot(1:24,SxV_N(1,:),'b',1:24,SxV_N(2,:),'r',1:24,SxV_N(3,:),'g');ylabel('%');xlabel('door#');
    subplot(3,2,6);plot(1:10,SxV_N(:,12),'b');ylim([0 50]);xlabel('% of goal visit');ylabel('segment');

end



