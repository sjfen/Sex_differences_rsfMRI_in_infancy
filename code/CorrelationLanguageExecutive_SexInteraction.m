%Spearman's rank-order correlation using functional
%connectivity and behavior datasets for between group (females vs males)
%comparison

clc
clear

%home directory
userdir = 'GITHUB_Directory';

%age string
ages_str0 = {'0YR', '1YR', '2YR'};
ages_str = {'neonate','oneyear','twoyear','fouryear', 'sixyear'};
ages_str2 = {'Neo', 'One', 'Two'};

%regions
region = {'ITGR', 'MTGL'};

for r = 1:length(region)
    % Set group
    for i = 1:3

        group=i;

        if group==1 % neonate
            aa = group;

            age_str0 = ages_str0{aa};
            age_str = ages_str{aa};
            age_str2 = ages_str2{aa};

            %load datafile
            data = readtable([userdir filesep 'data' filesep region{r} '_FC_Behavior_' age_str0],'Delimiter',',');

            BehaveNeo = data(:,end-8:end);

            %split by sex: behavior
            FemaleBehave_neonate = BehaveNeo(data.Sex == 1,:);
            MaleBehave_neonate = BehaveNeo(data.Sex == -1,:);

            if region{r} == 'ITGR'
                clustnum = [2 3];
            elseif region{r} == 'MTGL'
                clustnum = [2];
            end

            %residual extracted functional connectivity signal
            FC = table2array(data(:,clustnum));
            %split by sex: residual FC
            eval(['FemaleResidConnDat' region{r} '_' age_str ' = FC(data.Sex == 1,:)']);
            eval(['MaleResidConnDat' region{r} '_' age_str ' = FC(data.Sex == -1,:)']);

        elseif group==2 % oneyear
            aa = group;

            age_str0 = ages_str0{aa};
            age_str = ages_str{aa};
            age_str2 = ages_str2{aa};

            %load datafile
            data = readtable([userdir filesep 'data' filesep region{r} '_FC_Behavior_' age_str0],'Delimiter',',');

            BehaveOne = data(:,end-8:end);

            %split by sex: behavior
            FemaleBehave_oneyear = BehaveOne(data.Sex == 1,:);
            MaleBehave_oneyear = BehaveOne(data.Sex == -1,:);

            if region{r} == 'ITGR'
                clustnum = [2 3];
            elseif region{r} == 'MTGL'
                clustnum = [2 3];
            end

            %residual extracted functional connectivity signal
            FC = table2array(data(:,clustnum));
            %split by sex: residual FC
            eval(['FemaleResidConnDat' region{r} '_' age_str ' = FC(data.Sex == 1,:)']);
            eval(['MaleResidConnDat' region{r} '_' age_str ' = FC(data.Sex == -1,:)']);

        elseif group==3 %twoyear
            aa = group;

            age_str0 = ages_str0{aa};
            age_str = ages_str{aa};
            age_str2 = ages_str2{aa};

            %load datafile
            data = readtable([userdir filesep 'data' filesep region{r} '_FC_Behavior_' age_str0],'Delimiter',',');

            BehaveTwo = data(:,[end-8:end-5,end-1:end]); %not using 1YR language scores

            %split by sex: behavior
            FemaleBehave_twoyear= BehaveTwo(data.Sex == 1,:);
            MaleBehave_twoyear = BehaveTwo(data.Sex == -1,:);
            if region{r} == 'ITGR'
                clustnum = [2]; %make sure to zero out!
            elseif region{r} == 'MTGL'
                clustnum = [2 3];
            end

            %residual extracted functional connectivity signal
            FC = table2array(data(:,clustnum));
            %split by sex: residual FC
            eval(['FemaleResidConnDat' region{r} '_' age_str ' = FC(data.Sex == 1,:)']);
            eval(['MaleResidConnDat' region{r} '_' age_str ' = FC(data.Sex == -1,:)']);

        end

        for c =1:length(clustnum) %for each cluster used

            %FC correlate with behavior (Spearman's)
            %ANX (anxiety), BRIEFWM (working memory), BRIEFISCI (inhibitory
            %control), SBABIQ (intelligence), SBFIS (working memory), RL1
            %(receptive language 1yr), EL1 (expressive language 1yr), RL2
            %(receptive language 2yr), EL2 (expressive language 2yr)

            eval(['BhvAge = Behave' age_str2 ';'])
            for b = 1:size(BhvAge ,2)
                eval([ 'idxF = ((table2array(FemaleBehave_' age_str '(:,b)))~=0);']);
                eval([ 'idxM = ((table2array(MaleBehave_' age_str '(:,b)))~=0);']);
                eval(['idyF = (FemaleResidConnDat' region{r} '_' age_str '(:,c))~=0;']);
                eval(['idyM = (MaleResidConnDat' region{r} '_' age_str '(:,c))~=0;']);

                %interactions
                if sum(b==[1:5]) || sum(b==[6:7]) || sum(b==[8:9])

                    modelspec = 'Behavior ~ Sex + FC + Sex:FC';

                    eval(['FCtot=[FemaleResidConnDat' region{r} '_' age_str '([idxF & idyF],c); MaleResidConnDat' region{r} '_' age_str '([idxM & idyM],c)];']);
                    eval(['Sex=[ones(size(FemaleResidConnDat' region{r} '_' age_str '([idxF & idyF],c))); zeros(size(MaleResidConnDat' region{r} '_' age_str '([idxM & idyM],c)))];']);
                    eval(['Behavior=[table2array(FemaleBehave_' age_str '([idxF & idyF],b)); table2array(MaleBehave_' age_str '([idxM & idyM],b))];']);

                    tbl = table(FCtot,Sex,Behavior,'VariableNames', {'FC','Sex','Behavior'});

                    mdl = fitlm(tbl,modelspec);

                    eval(['interResid_pval_' region{r} '_' age_str '_clust' num2str(c) '_behav(b,:)= mdl.Coefficients.pValue(4);']);
                    eval(['interResid_est_' region{r} '_' age_str '_clust' num2str(c) '_behav(b,:)= mdl.Coefficients.Estimate(4);']);
                    eval(['interResid_tstat_' region{r} '_' age_str '_clust' num2str(c) '_behav(b,:)= mdl.Coefficients.tStat(4);']);

                end

            end

            eval(['table_FM_' region{r} age_str2 'Interact_clust' num2str(c) ' = table(interResid_est_' region{r} '_' age_str '_clust' num2str(c) '_behav, interResid_tstat_' region{r} '_' age_str '_clust' num2str(c) '_behav, interResid_pval_' region{r} '_' age_str '_clust' num2str(c) '_behav, ''VariableNames'', {''estimate'', ''tstat'', ''pval''});']);

        end
    end
end
