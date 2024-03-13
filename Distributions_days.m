%function [GxD_N,VxD_N,NxD_serial,D_NsegDist,MxDxTxG_PL,MxDxTxG_N,MxDxTxG_D,MxDxT_Nseg,MxD_NsegDist,MxD_door_visits,MxD_n_serial] = Distributions_days(MTSBELDG,MTBl)
%
%

function [GxD_N,VxD_N,NxD_serial,D_NsegDist,MxDxTxG_PL,MxDxTxG_N,MxDxTxG_D,MxDxT_Nseg,MxD_NsegDist,MxD_door_visits,MxD_n_serial] = Distributions_days(MTSBELDG,MTBl)

n_trials = 10;                                              %number of trials per day
n_days = length(unique(floor((MTBl(:,2)-1)/n_trials)));         %day1-->trial11:20, day2-->trial21:30
n_mice = length(unique(MTBl(:,1)));

%mattrices
MxDxTxG_PL = nan(n_mice,n_days,n_trials,25);
MxDxTxG_N = nan(n_mice,n_days,n_trials,25);
MxDxTxG_D = nan(n_mice,n_days,n_trials,25);
MxDxT_Nseg = nan(n_mice,n_days,n_trials);
MxD_NsegDist = nan(n_mice,n_days,50);
MxD_door_visits = zeros(n_mice,n_days,24);
MxD_n_serial = zeros(n_mice,n_days,20);
door_gaps = -12:1:12;
for ii = 1:n_days
    for mice_ii = 1:n_mice

        %number of segment by gap size, and path length by gap size
        for jj = 1:25
            for kk = 1:n_trials
                ndx = find(MTSBELDG(:,8)==door_gaps(jj) & MTSBELDG(:,2)==kk+ii*n_trials & MTSBELDG(:,1)==mice_ii);
                MxDxTxG_N(mice_ii,ii,kk,jj) = length(ndx); %per trial
                if ~isempty(ndx)
                    MxDxTxG_PL(mice_ii,ii,kk,jj) = mean(MTSBELDG(ndx,6));
                    MxDxTxG_D(mice_ii,ii,kk,jj) = mean(MTSBELDG(ndx,7));
                end                
            end
        end

        %number of segment per trial
        for jj = 1:n_trials
            ndx = find(MTSBELDG(:,2)==(ii*n_trials+jj) & MTSBELDG(:,1)==mice_ii);
            MxDxT_Nseg(mice_ii,ii,jj) = length(ndx);
        end
        MxD_NsegDist(mice_ii,ii,:) = histcounts(MxDxT_Nseg(mice_ii,ii,:),0.5:1:50.5);
        
        %door visits across days
        ndx = find(MTSBELDG(:,1)==mice_ii & MTSBELDG(:,2)>ii*n_trials & MTSBELDG(:,2)<=(ii+1)*n_trials);
        tmp = histcounts(MTSBELDG(ndx,5),-0.5:1:23.5);
        MxD_door_visits(mice_ii,ii,:) = [tmp(14:end) tmp(1:13)];

        %n_serial across days
        ndx = find(MTBl(:,1)==mice_ii & MTBl(:,2)>ii*n_trials & MTBl(:,2)<=(ii+1)*n_trials);
        MxD_n_serial(mice_ii,ii,:) = sum(MTBl(ndx,3:end),1);

    end
end


GxD_PL = squeeze(mean(MxDxTxG_PL,[1 3],'omitnan'))';
GxD_N = squeeze(mean(MxDxTxG_N,[1 3],'omitnan'))';
GxD_N = GxD_N./repmat(sum(GxD_N,1),size(GxD_N,1),1)*100;
VxD_N = squeeze(mean(MxD_door_visits,1,'omitnan'))';
VxD_N = VxD_N./repmat(sum(VxD_N,1),size(VxD_N,1),1)*100;
D_NsegDist = squeeze(mean(MxD_NsegDist,1))'/10*100;
NxD_serial = squeeze(mean(MxD_n_serial(:,:,1:20),1))';
NxD_serial = NxD_serial./repmat(sum(NxD_serial,1),size(NxD_serial,1),1)*100;

if nargout == 0
    D1 = 1;
    D2 = 6:15;
    figure;
    subplot(5,2,1);imagesc(GxD_PL);ylabel('seg length');xlabel('days');
    subplot(5,2,2);plot(-12:1:12,mean(GxD_PL(:,D1),2,'omitnan'),'b',-12:1:12,mean(GxD_PL(:,D2),2,'omitnan'),'r');ylabel('path length (cm)');xlabel('seg length');
    subplot(5,2,3);imagesc(GxD_N);ylabel('seg length');xlabel('days');
    subplot(5,2,4);plot(-12:1:12,mean(GxD_N(:,D1),2),'b',-12:1:12,mean(GxD_N(:,D2),2),'r');ylabel('% of segments');xlabel('seg length');
    subplot(5,2,5);imagesc(VxD_N);ylabel('door#');xlabel('days');
    subplot(5,2,6);plot(-11:12,mean(VxD_N(:,D1),2)/sum(mean(VxD_N(:,D1),2))*100,'b',-11:12,mean(VxD_N(:,D2),2)/sum(mean(VxD_N(:,D2),2))*100,'r');ylabel('% of visits');xlabel('door#');
    subplot(5,2,7);imagesc(NxD_serial);xlabel('days');ylabel('serial bout length');
    subplot(5,2,8);plot(1:20,mean(NxD_serial(:,D1),2),'b',1:20,mean(NxD_serial(:,D2),2),'r');xlabel('serial bout length');ylabel('% of serial bouts');
    subplot(5,2,9);imagesc(D_NsegDist);xlabel('days');ylabel('trial length');
    subplot(5,2,10);plot(1:50,mean(D_NsegDist(:,D1),2),'b',1:50,mean(D_NsegDist(:,D2),2),'r');xlabel('trial length (seg#)');ylabel('% of trials');
end



