clear, clc
tic

model_file = 'Cth_pyruvate_model';
mech_file = 'ctherm_pyruvate_mechanism';
data_file = 'Ctherm_pyruvate_data';



substrate = {'A','B','C','D','E','F','G','H','I','J','L','M'};
product = {'P','Q','R','S','T','U','V','W','X','Y','Z','N'};
inhibitor = {'a','b','c','d','e','f','g','h','p','j','l','m'};
activator = {'p','q','r','s','t','u','v','w','x','y','z'};

model = modcompile(model_file,mech_file,data_file);


MMparam = param_gen(model);

MM_final = struct('RXN',{''});
MM_symbolic = struct('RXN',{''});
for gg = 1:length(MMparam)
    ne = MMparam(gg).ne;
    np = MMparam(gg).np;
    param = MMparam(gg).param;

eq = cell(ne,ne);

for i = 1:ne
    for j = 1:ne
        if j==1
            e=j;
            for k = 1:np
                if k == 1
                    p = 2*k+(2*i-2);
                    if p > numel(param)
                        p = p - numel(param);
                    end
                    eq(1,i) = strcat(eq(1,i),param(p));
                    e = e + 1;
                    if e > ne
                        e = e - ne;
                    end
                else
                    p = 2*k+(2*i-2);
                    if p > numel(param)
                        p = p - numel(param);
                    end
                    eq(1,i) = strcat(eq(1,i),'*',param(p));
                    if e > ne
                        e = e-ne;
                    end
                end
            end
        elseif j == ne
            e = j;
            for k = 1:np
                if k == 1
                    p = 2*k+1+(2*i-2);
                    if p > numel(param)
                        p = p - numel(param);
                    end
                eq(j,i) = strcat(eq(j,i),param(p));
                e = e + 1;
                    if e > ne
                        e = e-ne;
                    end
                else
                    p = 2*k+1+(2*i-2);
                    if p > numel(param)
                        p = p - numel(param);
                    end
                    eq(ne,i) = strcat(eq(j,i),'*',param(p));
                    e = e + 1;
                    if e > ne
                        e = e-ne;
                    end
                end
            end
        else
            e = 1;
            for k = 1:(ne-j)
                if k == 1
                    p = 2*k+(2*i-2);
                    if p > numel(param)
                        p = p - numel(param);
                    end
                    eq(j,i) = strcat(eq(j,i),param(p));
                    e = e + 1;
                    if e > ne
                        e = e - ne;
                    end
                else
                    p = 2*k+(2*i-2);
                    if p > numel(param)
                        p = p - numel(param);
                    end
                    eq(j,i) = strcat(eq(j,i),'*',param(p));
                    e = e + 1;
                    if e > ne
                        e = e-ne;
                    end
                end
            end
            for k = e:np
                p = 2*k+1+(2*i-2);
                    if p > numel(param)
                        p = p - numel(param);
                    end
                eq(j,i) = strcat(eq(j,i),'*',param(p));
                e = e + 1;
                if e > ne
                    e = e - ne;
                end
            end
        end
        end
        

end


denom = [];
combine_denom = [];
kk = 1;
for i = 1:length(MMparam(gg).participant)
    if i == 1
        combine_denom = strcat(combine_denom,'isempty(strfind(eq{i,j}, ''',MMparam(gg).participant{i},'''))');
    else%if i ~= 1 && isempty(find(strcmp(MMparam(gg).participant{i},'p')))
        combine_denom = strcat(combine_denom,' && isempty(strfind(eq{i,j}, ''',MMparam(gg).participant{i},'''))');
    end
end

for i = 1:size(eq,1)
    for j = 1:size(eq,2)
        if eval(combine_denom) == 1
            if isempty(denom)
                denom = strcat(denom, eq{i,j});
            else
                denom = strcat(denom,'+',eq{i,j});
            end
        end
    end
end
        
denom = strsplit(denom,'+');
denom = denom';

if ~isempty(MMparam(gg).inhib_param)
    noinhibdim = size(eq,2);
    for i = 1:(length(MMparam(gg).inhib_param)/2)
        for j = 1:size(eq,1)
            for k = 1:size(eq,2)
                eq{j,k} = strcat(eq{j,k},'*',MMparam(gg).inhib_param{i*2});
            end 
        end
        for m = 1:size(denom,1)
            denom{m,1} = strcat(denom{m,1},'*',MMparam(gg).inhib_param{i*2});
        end
    end
        ll = 1;
        jj = 1;
        if ~isempty(model.kinetic(gg).c_in)
            for j = 1:length(model.kinetic(gg).c_in)
                eq{1,(noinhibdim+jj)} = [];
                for k = 1:size(eq,1)
                    eq{k,noinhibdim+jj} = strcat(eq{k,1},'*',MMparam(gg).inhib_param{ll});
                    split_inhib = strsplit(MMparam(gg).inhib_param{ll},'*');
                    split_val = sscanf(split_inhib{1},'k%d');
                    eq{k,noinhibdim+jj} = strrep(eq{k,noinhibdim+jj},strcat('*k',num2str(split_val+1)),'');
                end
                ll = ll+2;
                jj = jj + 1;
            end
        end
        if ~isempty(model.kinetic(gg).uc_in)
            
            for j = 1:length(model.kinetic(gg).uc_in)
                eq{1,(noinhibdim+jj)} = [];
                for k = 1:size(eq,1)
                    eq{k,(noinhibdim+jj)} = strcat(eq{k,2},'*',MMparam(gg).inhib_param{ll});
                    split_inhib = strsplit(MMparam(gg).inhib_param{ll},'*');
                    split_val = sscanf(split_inhib{1},'k%d');
                    eq{k,noinhibdim+jj} = strrep(eq{k,noinhibdim+jj},strcat('*k',num2str(split_val+1)),'');
                end
                ll = ll + 2;
                jj = jj + 1;
            end
        end
        if ~isempty(model.kinetic(gg).nc_in)            
            for j = 1:length(model.kinetic(gg).nc_in)
                eq{1,(noinhibdim+jj)} = [];
                eq{1,(noinhibdim+jj+1)} = [];                
                for k = 1:size(eq,1)
                    eq{k,(noinhibdim+jj)} = strcat(eq{k,1},'*',MMparam(gg).inhib_param{ll});
                    split_inhib = strsplit(MMparam(gg).inhib_param{ll},'*');
                    split_val = sscanf(split_inhib{1},'k%d');
                    eq{k,noinhibdim+jj} = strrep(eq{k,noinhibdim+jj},strcat('*k',num2str(split_val+1)),'');                    
                    eq{k,(noinhibdim+jj+1)} = strcat(eq{k,2},'*',MMparam(gg).inhib_param{ll+2});
                    split_inhib = strsplit(MMparam(gg).inhib_param{ll+2},'*');
                    split_val = sscanf(split_inhib{1},'k%d');
                    eq{k,noinhibdim+jj+1} = strrep(eq{k,noinhibdim+jj+1},strcat('*k',num2str(split_val+1)),'');                    
                end
                ll = ll + 4;
                jj = jj + 2;                
            end
        end    

end
if ~isempty(MMparam(gg).act_param)
    noactdim = size(eq,2);
    for i = 1:(length(MMparam(gg).act_param)/2)
        for j = 1:size(eq,1)
            for k = 1:size(eq,2)
                eq{j,k} = strcat(eq{j,k},'*',MMparam(gg).act_param{i*2});
            end 
        end
        for m = 1:size(denom,1)
            denom{m,1} = strcat(denom{m,1},'*',MMparam(gg).act_param{i*2});
        end
    end
        ll = 1;
        jj = 1;
        if ~isempty(model.kinetic(gg).act)
            for j = 1:length(model.kinetic(gg).act)
                eq{1,(noactdim+jj)} = [];
                for k = 1:size(eq,1)
                    eq{k,noactdim+jj} = strcat(eq{k,1},'*',MMparam(gg).act_param{ll});
                    split_act = strsplit(MMparam(gg).act_param{ll},'*');
                    split_val = sscanf(split_act{1},'k%d');
                    eq{k,noactdim+jj} = strrep(eq{k,noactdim+jj},strcat('*k',num2str(split_val+1)),'');
                    eq{k,noactdim+jj} = strrep(eq{k,noactdim+jj},'p*','');
                end
                ll = ll+2;
                jj = jj + 1;
            end
        end    
end
    

Mparam(gg).eq = eq;


for i = 1:size(eq,1)
    if i == 1
        kcat1 = strcat(param(numel(param)-1),'*',eq{i,ne});
    else
        kcat1 = strcat(kcat1,'+',param(numel(param)-1),'*',eq{i,ne});
    end
end
for i = 1:size(eq,1)
    if i == 1
        kcat2 = strcat(eq{i,1},'*',param(numel(param)));
    else
        kcat2 = strcat(kcat2,'+',eq{i,1},'*',param(numel(param)));
    end
end

aaa=cell(1,1);
aaa{1,1} = kcat1;
bbb=cell(1,1);
bbb{1,1}=kcat2;
A=str2sym(string(aaa));
B=str2sym(string(bbb));
syms RR
kcatterms = solve(RR == A - B,RR);
kcatterms = char(kcatterms);
kcatterms = strsplit(kcatterms,' - ');
kcatterms = kcatterms';
MMparam(gg).kcatterms = kcatterms;

% for i = 1:size(eq,1)
%     for j = 1:size(eq,2)
%         if i == 1 && j==1
%             denom = eq{i,j};
%         else
%             denom = strcat(denom,'+',eq{i,j});
%         end
%     end
% end

denom2 = [];
MM=[];
%denom = [];

combine_denom = [];
kk = 1;
for i = 1:length(MMparam(gg).participant)
    if i == 1
        combine_denom = strcat(combine_denom,'isempty(strfind(eq{i,j}, ''',MMparam(gg).participant{i},'''))');
    elseif i ~= 1 && isempty(find(strcmp(MMparam(gg).participant{i},'p')))
        combine_denom = strcat(combine_denom,' && isempty(strfind(eq{i,j}, ''',MMparam(gg).participant{i},'''))');
    end
end

for i = 1:size(eq,1)
    for j = 1:size(eq,2)
        if ~isempty(find(strcmp(eq{i,j},denom)))
            holder_var = 1;
%             if isempty(denom)
%                 denom = strcat(denom,eq{i,j});
%             else
%                 denom = strcat(denom,'+',eq{i,j});
%             end
        else
            if isempty(MM)
                MM = strcat(MM,eq{i,j});
            else
                MM = strcat(MM,'+',eq{i,j});
            end
        end
    end
end
for imb = 1:length(denom)
    if imb == 1
        denom2 = denom{imb};
    else
        denom2 = strcat(denom2,'+',denom{imb});
    end
end
MMparam(gg).denom = denom2;
MMparam(gg).MM = MM; 


a = [];
enu = [];
participant = MMparam(gg).participant;
for i = 1:length(participant)
    a = strcat(a,participant{i});
end

for i = 1:length(participant)
    enu{i,1} = combnk(a,i);
end
MMparam(gg).enu = enu;

separate = strsplit(MM,'+');
clear len_enu
for i = 1:length(enu)
    len_enu(i) = size(char(enu(i)),1);
end

MM_formula = cell(sum(len_enu),2);
d = 1;

for i = 1:length(enu)
    current_enu = char(enu(i));
    size1_enu = char(enu(1));
    
    for h = 1:size(current_enu,1)
        combine_enu = [];
        for j = 1:size(current_enu,2)
            if j == 1
                combine_enu = strcat(combine_enu,'~isempty(strfind(separate{kk}, ''',current_enu(h,j),'''))');
            else
                combine_enu = strcat(combine_enu,' && ~isempty(strfind(separate{kk}, ''',current_enu(h,j),'''))');
            end
        end
        for ll = 1:size(size1_enu,1);
            leave_out = strfind(current_enu(h,:),size1_enu(ll));
            if isempty(leave_out)
                combine_enu = strcat(combine_enu,' && isempty(strfind(separate{kk}, ''',size1_enu(ll),'''))');
            end
        end        
        MM_formula{d,1} = current_enu(h,:);

            for kk = 1:length(separate)
                if eval(combine_enu) == 1
                    if isempty(MM_formula{d,2})
                        MM_formula{d,2} = strcat(MM_formula{d,2},separate(kk));
                    else
                        MM_formula{d,2} = strcat(MM_formula{d,2},'+',separate(kk));
                    end
                end
            end
            d = d + 1;
        end
end
put_form = size(MM_formula,1);
if ~isempty(find(strcmp('p',participant)))
    
    MM_formula{put_form+1,1} = '';
    pull_param = strsplit(MMparam(gg).MM,'+');
    grg = 1;
    for imdb = 1:length(pull_param)
        end_g = 0;
        for gbg = 1:length(substrate)
            if ~isempty(strfind(pull_param{imdb},substrate{gbg}))
                end_g = end_g + 1;
            end
        end
        for gbg = 1:length(product)
            if ~isempty(strfind(pull_param{imdb},product{gbg}))
                end_g = end_g + 1;
            end
        end
        for gbg = 1:length(inhibitor)
            if ~isempty(strfind(pull_param{imdb},inhibitor{gbg}))
                end_g = end_g + 1;
            end
        end
        for gbg = 1:length(activator)
            if ~isempty(strfind(pull_param{imdb},activator{gbg}))
                end_g = end_g + 1;
            end
        end
        if end_g == 0
            if grg == 1;
                MM_formula{put_form+1,2} = pull_param{imdb};
                grg = grg + 1;
            else
                MM_formula{put_form+1,2} = strcat(MM_formula{put_form+1,2},'+',pull_param{imdb});
            end
        end
    end
end
            
                
    

if numel(participant)>1
    

tot_permu = size(MM_formula,1);
    for i = 1:tot_permu
        if isempty(MM_formula{(tot_permu-i+1),2})
            MM_formula((tot_permu-i+1),:) = [];
        end
    end
    %MM_final{gg} = MM_formula;
    MM_final(gg).RXN = MMparam(gg).id;
    MM_final(gg).kcatterms = MMparam(gg).kcatterms;
    MM_final(gg).MM = MM_formula;
    MM_final(gg).denom = MMparam(gg).denom;
else
    MMparam(gg).denom = 'k2 + k3';
    MMparam(gg).MM = cell(1,2);
    MMparam(gg).MM{1,1} = cellstr('A');
    MMparam(gg).MM{1,2} = cellstr('k1*A');
    MMparam(gg).kcatterms = cellstr('k3');
    MM_formula = cell(1,1);
    MM_formula{1} = 'k1*A';
    MM_final(gg).RXN = MMparam(gg).id;
    MM_final(gg).kcatterms = MMparam(gg).kcatterms;
    MM_final(gg).MM = MMparam(gg).MM;
    MM_final(gg).denom = MMparam(gg).denom;
end
    
    
    
    MM_symbolic(gg).RXN = MM_final(gg).RXN;
    for i = 1:size(MM_final(gg).MM,1)
        MM_symbolic(gg).MM{i,1} = MM_final(gg).MM{i,1};
        MM_symbolic(gg).MM{i,2} = str2sym(string(MM_final(gg).MM{i,2}));
    end
    MM_symbolic(gg).denom = str2sym(MM_final(gg).denom);
    for i = 1:length(MM_final(gg).kcatterms)
        MM_symbolic(gg).kcat{i,1} = str2sym(MM_final(gg).kcatterms{i});
    end
    
    if length(MM_symbolic) == 37
        cchchc = 1;
    end
    for i = 1:size(MM_symbolic(gg).MM,1)
        MM_symbolic(gg).MM{i,2} = MM_symbolic(gg).denom/MM_symbolic(gg).MM{i,2};
    end
    for i = 1:length(MM_symbolic(gg).kcat)
        MM_symbolic(gg).kcatterms{i,1} = MM_symbolic(gg).kcat{i,1}/MM_symbolic(gg).denom;
    end
end

%%%%%%%%%%%%
for i = 1:length(MM_symbolic)
    ll = 1;
    bb = 1;
    clear inc_sub
    clear inc_prod
        for k = 1:length(MMparam(i).participant)
            if ~isempty(find(strcmp(MMparam(i).participant{k},substrate)))
                inc_sub(ll,1) = MMparam(i).participant{k};
                ll = ll + 1;
            end
            if ~isempty(find(strcmp(MMparam(i).participant{k},product)))
                inc_prod(bb,1) = MMparam(i).participant{k};
                bb = bb + 1;
            end
        end
        if exist('inc_prod') == 0
            inc_prod = [];
        end
        if exist('inc_sub') == 0
            inc_sub = [];
        end
        for k = 1:size(MM_symbolic(i).MM,1)
            for l = 1:length(substrate)
                subcomp(l) = isempty(strfind(MM_symbolic(i).MM{k,1},substrate{l}));
            end
            for l = 1:length(product)
                prodcomp(l) = isempty(strfind(MM_symbolic(i).MM{k,1},product{l}));
            end
            for l = 1:length(inhibitor)
                inhibcomp(l) = isempty(strfind(MM_symbolic(i).MM{k,1},inhibitor{l}));
            end
            for l = 1:length(activator)
                actcomp(l) = isempty(strfind(MM_symbolic(i).MM{k,1},activator{l}));
            end
            subsum = 0;
            prodsum = 0;
            for l = 1:length(inc_sub)
                if ~isempty(strfind(inc_sub(l),MM_symbolic(i).MM{k,1}))
                    subsum = subsum+1;
                end
            end
            for l = 1:length(inc_prod)
                if ~isempty(strfind(inc_prod(l),MM_symbolic(i).MM{k,1}))
                    prodsum = prodsum+1;
                end
            end
            if subsum == length(inc_sub) && sum(prodcomp) == 12 && sum(inhibcomp) == 12 && sum(actcomp) == 11
                MM_symbolic(i).kcatterms{1} = MM_symbolic(i).kcatterms{1}*MM_symbolic(i).MM{k,2};
            end
            if prodsum == length(inc_prod) && sum(subcomp) == 12 && sum(inhibcomp) == 12 && sum(actcomp) == 11
                MM_symbolic(i).kcatterms{2} = MM_symbolic(i).kcatterms{2}*MM_symbolic(i).MM{k,2};
            end
        end
end


%%%%%%%%%%%%


%load MM_symbolic_to47.mat

for i = 1:size(MM_symbolic,2)
%{    
    if i == 38| i == 47
        nb = 10;
    elseif i ==49
        nb = 3;
    elseif i == 50 | i == 52 | i == 54| i== 55 | i == 56 | i == 57 |i == 46
        nb = 5;
    else
        nb = size(MM_symbolic(i).MM,1);
end
%}
    for j = 1:size(MM_symbolic(i).MM,1)
        
%         if i == 38 && j == 18
%             aba = strrep(char(MM_symbolic(i).MM{j,2}),'*1i)','');
%             aba = strrep(aba,'-(','');
%             MM_symbolic(i).MM{j,2} = sym(aba);
%         end
%         if i == 38 && j >= 28
%             aba = strrep(char(MM_symbolic(i).MM{j,2}),'*1i)','');
%             aba = strrep(aba,'-(','');
%             aba = strrep(aba,'exp(-1)*','');
%             MM_symbolic(i).MM{j,2} = sym(aba);
%         end
        
        working_mm = char(MM_symbolic(i).MM{j,2});
        working_mm = strsplit(working_mm,'/');
        working_num = char(working_mm(1));
        working_num = strrep(working_num,'(','');
        working_num = strrep(working_num,')','');
        working_num = strsplit(working_num,' + ');
        working_denom = char(working_mm(2));
        working_denom = strrep(working_denom,'(','');
        working_denom = strrep(working_denom,')','');
        working_denom = strsplit(working_denom,'+');
        for_check = char(working_num(1));
        check_num = strsplit(for_check,'*');
        for k = 1:length(check_num)
            if strcmp(check_num(k),'k1') ~= 1 && strcmp(check_num(k),'k2') ~= 1 && strcmp(check_num(k),'k3') ~= 1
                    store = [];
                    store2 = [];
                    for m = 1:length(working_num)
                    num_eval = strfind(working_num(m),check_num(k));
                    if ~isempty(num_eval{1})
                        store(m) = 1;
                    else
                        store(m) = 0;
                    end
                    end
                    for m = 1:length(working_denom)
                        denom_eval = strfind(working_denom(m),check_num(k));
                        if ~isempty(denom_eval{1})
                            store2(m) = 1;
                        else
                            store2(m) = 0;
                        end
                    end
                    if sum(store) == length(working_num) && sum(store2) == length(working_denom)
                        working_mm = strrep(working_mm,strcat('*',check_num(k)),'');
                    end
            end
                

        end
       %working_mm = strcat(working_mm(1),'/',working_mm(2));
       
       inhib_mm = working_mm;%strsplit(working_mm,'/');
       n = 2;
       for m = 1:length(MMparam(i).inhib_param)/2
           inhib_check = MMparam(i).inhib_param(n);
           inhib_f = char(MMparam(i).inhib_param(n-1));
           emptcheck = strfind(inhib_mm(1),inhib_check{1}); 
           if ~isempty(emptcheck{1,1})
               inhib_mm(1) = strrep(inhib_mm(1),strcat('*',inhib_check),'');
               inhib_f = strsplit(inhib_f,'*');
               inhib_mm(2) = strrep(inhib_mm(2),strcat('*',inhib_f(1)),'');
               inhib_mm(1) = strcat('kI(',num2str(m),'*',inhib_mm(1));
           end
           n = n + 2;
       end
       inhib_mm = strcat(inhib_mm(1),'/',inhib_mm(2));
       inhib_mm = insertAfter(inhib_mm,'k','(');
       inhib_mm = insertBefore(inhib_mm,')',')');
       inhib_mm = insertBefore(inhib_mm,' +',')');
       inhib_mm = insertBefore(inhib_mm,'*',')');
       substrate = {'A','B','C','D','E','F','G','H','I','J','L','M'};
       product = {'P','Q','R','S','T','U','V','W','X','Y','Z','N'};
       inhibitor = {'a','b','c','d','e','f','g','h','i','j','l','m'};
       activator = {'p','q','r','s','t','u','v','w','x','y','z'};

       for mm = 1:length(substrate)
           inhib_mm = strrep(inhib_mm,strcat(substrate(mm),')'),substrate(mm));
       end
       for mm = 1:length(product)
           inhib_mm = strrep(inhib_mm,strcat(product(mm),')'),product(mm));
       end
       for mm = 1:length(inhibitor)
           inhib_mm = strrep(inhib_mm,strcat(inhibitor(mm),')'),inhibitor(mm));
       end
       inhib_mm = strrep(inhib_mm,'k(I','kI');
       
       act_mm = inhib_mm;%strsplit(working_mm,'/');
       n = 2;
       if ~isempty(MMparam(i).act_param)
           if ~isempty(strfind(act_mm,MMparam(i).act_param{2}))
               act_mm = strcat('(kA(1)*',act_mm);
               act_sub = strrep(MMparam(i).act_param{2},'k','k(');
               act_sub = strrep(act_sub,'*',')*');
               act_mm = strrep(act_mm,act_sub,'p');
               act_mm = strrep(act_mm,'p)','p');
               act_mm = strrep(act_mm,'/',')/');
               act_sub2 = strrep(MMparam(i).act_param{1},'k','k(');
               act_sub2 = strcat(act_sub2,')');
               act_mm = strrep(act_mm,strcat('*',act_sub2),'');
               inhib_mm = act_mm;
           end
       end
%        act_mm = strcat(act_mm(1),'/',act_mm(2));
%        act_mm = insertAfter(act_mm,'k','(');
%        act_mm = insertBefore(act_mm,')',')');
%        act_mm = insertBefore(act_mm,' +',')');
%        act_mm = insertBefore(act_mm,'*',')');
%        substrate = {'A','B','C','D','E','F','G','H','I','J','L','M'};
%        product = {'P','Q','R','S','T','U','V','W','X','Y','Z','N'};
%        inhibitor = {'a','b','c','d','e','f','g','h','i','j','l','m'};
%        activator = {'p','q','r','s','t','u','v','w','x','y','z'};
% 
%        for mm = 1:length(substrate)
%            act_mm = strrep(act_mm,strcat(substrate(mm),')'),substrate(mm));
%        end
%        for mm = 1:length(product)
%            act_mm = strrep(act_mm,strcat(product(mm),')'),product(mm));
%        end
%        for mm = 1:length(inhibitor)
%            act_mm = strrep(act_mm,strcat(inhibitor(mm),')'),inhibitor(mm));
%        end
%        for mm = 1:length(activator)
%            act_mm = strrep(act_mm,strcat(activator(mm),')'),activator(mm));
%        end
%        act_mm = strrep(act_mm,'k(A','kA')
       
%{
       if i == 49 | i == 52
           inhib_mm = strrep(inhib_mm,'exp(-1)))*','');
           inhib_mm = strrep(inhib_mm,'))))/',')))/');
       end
       if i == 57
           inhib_mm = strrep(inhib_mm,'-','');
           inhib_mm = strrep(inhib_mm,'*1i)','');
       end
%}
       
       
       MM_symbolic(i).MM{j,2} = str2sym(string(inhib_mm));
       
    end
end


%%%%%%%


for i = 1:size(MM_symbolic,2)
    for j = 1:size(MM_symbolic(i).kcatterms,1)
        working_kcat = char(MM_symbolic(i).kcatterms{j,1});
        working_kcat = strsplit(working_kcat,'/');
        working_num = char(working_kcat(1));
        working_num = strrep(working_num,'(','');
        working_num = strrep(working_num,')','');
        working_num = strsplit(working_num,' + ');
        working_denom = char(working_kcat(2));
        working_denom = strrep(working_denom,'(','');
        working_denom = strrep(working_denom,')','');
        working_denom = strsplit(working_denom,'+');
        for_check = char(working_num(1));
        check_num = strsplit(for_check,'*');
        for k = 1:length(check_num)
            if strcmp(check_num(k),'k1') ~= 1 && strcmp(check_num(k),'k2') ~= 1 && strcmp(check_num(k),'k3') ~= 1
                    store = [];
                    store2 = [];
                    for m = 1:length(working_num)
                    num_eval = strfind(working_num(m),check_num(k));
                    if ~isempty(num_eval{1})
                        store(m) = 1;
                    else
                        store(m) = 0;
                    end
                    end
                    for m = 1:length(working_denom)
                        denom_eval = strfind(working_denom(m),check_num(k));
                        if ~isempty(denom_eval{1})
                            store2(m) = 1;
                        else
                            store2(m) = 0;
                        end
                    end
                    if sum(store) == length(working_num) && sum(store2) == length(working_denom)
                        working_kcat = strrep(working_kcat,strcat('*',check_num(k)),'');
                    end
            end
        end

        working_kcat = strcat(working_kcat(1),')/',working_kcat(2));
        working_kcat = insertAfter(working_kcat,'k','(');
        working_kcat = insertBefore(working_kcat,' +',')');
        working_kcat = insertBefore(working_kcat,'*',')');

        for mm = 1:length(substrate)
            working_kcat = strrep(working_kcat,strcat(substrate(mm),')'),substrate(mm));
            working_kcat = strrep(working_kcat,strcat(product(mm),')'),product(mm));
            working_kcat = strrep(working_kcat,strcat(inhibitor(mm),')'),inhibitor(mm));
            if mm < 12
                working_kcat = strrep(working_kcat,strcat(activator(mm),')'),activator(mm));
            end
        end
        working_kcat = strcat(working_kcat,')');
        %if i == 38 | i == 47 | i == 49 | i == 50 | i == 52 | i == 54 | i == 55 | i == 46 | i == 56 | i == 57
        %    working_kcat = strrep(working_kcat,'*exp(1))','');
        %    working_kcat = strrep(working_kcat,'*1i','');
        %end
        %if i == 46
        %   working_kcat = strrep(working_kcat,'*exp(1)','');
        %   working_kcat = strrep(working_kcat,'*1i',''); 
        %end
        MM_symbolic(i).kcatterms{j,1} = str2sym(string(working_kcat));
    end
end


for i = 1:size(model.ensemble.ne,1)
    dkdc{i,1} = model.kinetic(i).dkdc;
end
p = 1;
for i = 1:size(dkdc,1)
    currentrxn = dkdc{i};
    for j = 1:size(currentrxn,2)
        a = find(currentrxn(:,j));
        if ~isempty(a)
            kinetic_name{p,1} = model.metprop(a).metid;
            p = p + 1;
        else
            kinetic_name{p,1} = model.kinetic(i).id;
            p = p + 1;
        end
    end   
end
p = 1;
for i = 1:size(dkdc,1)
    for j = 1:model.kinetic(i).nk
        RXN_NAME{p,1} = model.kinetic(i).id;
        p = p + 1;
    end
end

%[~,metabID,~] = xlsread('metabolomics_ranges.xlsx','range','C2:C46');
%[concLB,~,~] = xlsread('metabolomics_ranges.xlsx','range','D2:D46');
%[concUB,~,~] = xlsread('metabolomics_ranges.xlsx','range','E2:E46');



kLB = xlsread('pyruvate_param_052521.xlsx','Sheet1','K2:K79');
kUB = xlsread('pyruvate_param_052521.xlsx','Sheet1','L2:L79');

for i = 1:length(model.kinetic)
    RxnName{i,1} = model.kinetic(i).id;
end
model.RxnName = RxnName;

%kLB = xlsread('Ali_BOUNDS_11_30_18.xlsx','BOUNDS','A2:A309');
%kUB = xlsread('Ali_BOUNDS_11_30_18.xlsx','BOUNDS','B2:B309');

exclude = cell(0,5);
solo = 1;
for i = 1:length(model.RxnName)
    kI = [];
    k = [];
    kA = [];
    vn = 1;
    param_1 = find(strcmp(model.RxnName{i},RXN_NAME));
    for j = 1:length(MMparam(i).param)
        MM_symbolic(i).kLB(j,1) = kLB(param_1(j));
        MM_symbolic(i).kUB(j,1) = kUB(param_1(j));
        MM_symbolic(i).metab_k{j,1} = kinetic_name{param_1(j),1};
        %if ~isempty(find(strcmp(kinetic_name{param_1(j),1},metabID)))
        %    MM_symbolic(i).kLB(j,1) = MM_symbolic(i).k(j,1)/(1000*concLB(find(strcmp(kinetic_name{param_1(j),1},metabID))));
        %    MM_symbolic(i).kUB(j,1) = MM_symbolic(i).k(j,1)/(1000*concUB(find(strcmp(kinetic_name{param_1(j),1},metabID))));
        %else
            %MM_symbolic(i).kLB(j,1) = MM_symbolic(i).k(j,1);
            %MM_symbolic(i).kUB(j,1) = MM_symbolic(i).k(j,1);
            if isempty(find(strcmp(kinetic_name{param_1(j),1},model.RxnName)))
                exclude{solo,1} = model.RxnName{i};
                exclude{solo,2} = 'substrate';
                exclude{solo,3} = j;
                exclude{solo,4} = MM_symbolic(i).metab_k{j,1};
                exclude{solo,5} = MMparam(i).param{j,1};
                solo = solo + 1;
            end
                
        %end
    end
    for j = 1:length(model.kinetic(i).c_in)
        MM_symbolic(i).kILB(vn,1) = kLB(param_1(length(MMparam(i).param)+vn));
        MM_symbolic(i).kIUB(vn,1) = kUB(param_1(length(MMparam(i).param)+vn));
        MM_symbolic(i).metab_kI{vn,1} = kinetic_name{param_1(length(MMparam(i).param)+vn),1};
        %if ~isempty(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+vn),1},metabID)))
        %    MM_symbolic(i).kILB(vn,1) = MM_symbolic(i).kI(vn,1)*(1000*concLB(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+vn),1},metabID))));
        %    MM_symbolic(i).kIUB(vn,1) = MM_symbolic(i).kI(vn,1)*(1000*concUB(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+vn),1},metabID))));
        %else
            %MM_symbolic(i).kILB(vn,1) = MM_symbolic(i).kI(vn,1)
            %MM_symbolic(i).kIUB(vn,1) = MM_symbolic(i).kI(vn,1)
            if isempty(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+vn),1},model.RxnName)))
                exclude{solo,1} = model.RxnName{i};
                exclude{solo,2} = 'inhibitor';
                exclude{solo,3} = vn;
                exclude{solo,4} = MM_symbolic(i).metab_kI{vn,1};
                exclude{solo,5} = MMparam(i).inhib_param{vn,1};
                solo = solo + 1;
            end
        %end
        vn = vn+1;
    end
    for j = 1:length(model.kinetic(i).uc_in)
        MM_symbolic(i).kILB(vn,1) = kLB(param_1(length(MMparam(i).param)+vn));
        MM_symbolic(i).kIUB(vn,1) = kUB(param_1(length(MMparam(i).param)+vn));
        MM_symbolic(i).metab_kI{vn,1} = kinetic_name{param_1(length(MMparam(i).param)+vn),1};
        %if ~isempty(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+vn),1},metabID)))
        %    MM_symbolic(i).kILB(vn,1) = MM_symbolic(i).kI(vn,1)*(1000*concLB(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+vn),1},metabID))));
        %    MM_symbolic(i).kIUB(vn,1) = MM_symbolic(i).kI(vn,1)*(1000*concUB(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+vn),1},metabID))));
        %else
            %MM_symbolic(i).kILB(vn,1) = MM_symbolic(i).kI(vn,1)
            %MM_symbolic(i).kIUB(vn,1) = MM_symbolic(i).kI(vn,1)
            if isempty(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+vn),1},model.RxnName)))
                exclude{solo,1} = model.RxnName{i};
                exclude{solo,2} = 'inhibitor';
                exclude{solo,3} = vn;
                exclude{solo,4} = MM_symbolic(i).metab_kI{vn,1};
                exclude{solo,5} = MMparam(i).inhib_param{vn,1};
                solo = solo + 1;
            end
        %end
        vn = vn+1;
    end
    kn = vn;
    for j = 1:length(model.kinetic(i).nc_in)
        MM_symbolic(i).kILB(kn,1) = kLB(param_1(length(MMparam(i).param)+vn));
        MM_symbolic(i).kIUB(kn,1) = kUB(param_1(length(MMparam(i).param)+vn));
        MM_symbolic(i).metab_kI{kn,1} = kinetic_name{param_1(length(MMparam(i).param)+vn),1};
        %if ~isempty(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+vn),1},metabID)))
        %    MM_symbolic(i).kILB(kn,1) = MM_symbolic(i).kI(kn,1)*(1000*concLB(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+vn),1},metabID))));
        %    MM_symbolic(i).kIUB(kn,1) = MM_symbolic(i).kI(kn,1)*(1000*concUB(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+vn),1},metabID))));
        %else
            %MM_symbolic(i).kILB(kn,1) = MM_symbolic(i).kI(kn,1)
            %MM_symbolic(i).kIUB(kn,1) = MM_symbolic(i).kI(kn,1)
            if isempty(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+vn),1},model.RxnName)))
                exclude{solo,1} = model.RxnName{i};
                exclude{solo,2} = 'inhibitor';
                exclude{solo,3} = j;
                exclude{solo,4} = MM_symbolic(i).metab_kI{vn,1};
                exclude{solo,5} = MMparam(i).inhib_param{vn,1};
                solo = solo + 1;
            end
        %end
        MM_symbolic(i).kILB(kn+1,1) = kLB(param_1(length(MMparam(i).param)+length(model.kinetic(i).nc_in)+vn));
        MM_symbolic(i).kIUB(kn+1,1) = kUB(param_1(length(MMparam(i).param)+length(model.kinetic(i).nc_in)+vn));
        MM_symbolic(i).metab_kI{kn+1,1} = kinetic_name{param_1(length(MMparam(i).param)+length(model.kinetic(i).nc_in)+vn),1};
        %if ~isempty(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+vn),1},metabID)))
        %    MM_symbolic(i).kILB(kn+1,1) = MM_symbolic(i).kI(kn+1,1)*(1000*concLB(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+length(model.kinetic(i).nc_in)+vn),1},metabID))));
        %    MM_symbolic(i).kIUB(kn+1,1) = MM_symbolic(i).kI(kn+1,1)*(1000*concUB(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+length(model.kinetic(i).nc_in)+vn),1},metabID))));
        %else
            %MM_symbolic(i).kILB(kn+1,1) = MM_symbolic(i).kI(kn+1,1)
            %MM_symbolic(i).kIUB(kn+1,1) = MM_symbolic(i).kI(kn+1,1)
            if isempty(find(strcmp(kinetic_name{param_1(length(MMparam(i).param)+vn),1},model.RxnName)))
                exclude{solo,1} = model.RxnName{i};
                exclude{solo,2} = 'inhibitor';
                exclude{solo,3} = j;
                exclude{solo,4} = MM_symbolic(i).metab_kI{vn,1};
                exclude{solo,5} = MMparam(i).inhib_param{vn,1};
                solo = solo + 1;
            end
        %end
        vn = vn+1;
        kn = kn + 2;
    end
end

if isfield(MM_symbolic,'kILB') == 0
    for i = 1:size(MM_symbolic,2)
        MM_symbolic(i).kILB = [];
    end
end
if isfield(MM_symbolic,'kIUB') == 0
    for i = 1:size(MM_symbolic,2)
        MM_symbolic(i).kIUB = [];
    end
end


for i = 1:size(MM_symbolic,2)
    
    %{
    if i == 38 | i == 47
        nb = 10;
    elseif i == 46
        nb = 5;
    elseif i == 49
        nb = 3
    else
        nb = size(MM_symbolic(i).MM,1);
    end
    %}
    for j = 1:size(MM_symbolic(i).MM,1)
        MMparameters_final(i).Km{j,1} = MM_symbolic(i).MM{j,1};
        for k = 1:length(substrate)
            if k == 1
                MMparameters_final(i).Km{j,2} = strrep(char(MM_symbolic(i).MM{j,2}),strcat(substrate{k},'*'),'');
                MMparameters_final(i).Km{j,2} = strrep(char(MMparameters_final(i).Km{j,2}),strcat('*',substrate{k}),'');
            else
                MMparameters_final(i).Km{j,2} = strrep(char(MMparameters_final(i).Km{j,2}),strcat(substrate{k},'*'),'');
                MMparameters_final(i).Km{j,2} = strrep(char(MMparameters_final(i).Km{j,2}),strcat('*',substrate{k}),'');
            end
        end
        for k = 1:length(product)
            if k == 1
                MMparameters_final(i).Km{j,2} = strrep(char(MMparameters_final(i).Km{j,2}),strcat(product{k},'*'),'');
                MMparameters_final(i).Km{j,2} = strrep(char(MMparameters_final(i).Km{j,2}),strcat('*',product{k}),'');
            else
                MMparameters_final(i).Km{j,2} = strrep(char(MMparameters_final(i).Km{j,2}),strcat(product{k},'*'),'');
                MMparameters_final(i).Km{j,2} = strrep(char(MMparameters_final(i).Km{j,2}),strcat('*',product{k}),'');
            end  
        end
        for k = 1:length(inhibitor)
            if k == 1
                MMparameters_final(i).Km{j,2} = strrep(char(MMparameters_final(i).Km{j,2}),strcat(inhibitor{k},'*'),'');
                MMparameters_final(i).Km{j,2} = strrep(char(MMparameters_final(i).Km{j,2}),strcat('*',inhibitor{k}),'');
            else
                MMparameters_final(i).Km{j,2} = strrep(char(MMparameters_final(i).Km{j,2}),strcat(inhibitor{k},'*'),'');
                MMparameters_final(i).Km{j,2} = strrep(char(MMparameters_final(i).Km{j,2}),strcat('*',inhibitor{k}),'');
            end
        end
        for k = 1:length(activator)
            if k == 1
                MMparameters_final(i).Km{j,2} = strrep(char(MMparameters_final(i).Km{j,2}),strcat(activator{k},'*'),'');
                MMparameters_final(i).Km{j,2} = strrep(char(MMparameters_final(i).Km{j,2}),strcat('*',activator{k}),'');
            else
                MMparameters_final(i).Km{j,2} = strrep(char(MMparameters_final(i).Km{j,2}),strcat(activator{k},'*'),'');
                MMparameters_final(i).Km{j,2} = strrep(char(MMparameters_final(i).Km{j,2}),strcat('*',activator{k}),'');
            end
        end
        klo = [MM_symbolic(i).kLB;MM_symbolic(i).kILB];
        kup = [MM_symbolic(i).kUB;MM_symbolic(i).kIUB];
        kinit = (klo+kup)/2;
        zzz = length(MMparam(i).param);
        Km_store = MMparameters_final(i).Km{j,2};
        for k = 1:length(MM_symbolic(i).kIUB)
            kk = strcat('kI(',num2str(k),')');
            if strfind(Km_store,kk)
                Km_store = strrep(Km_store,kk,strcat('k(',num2str(zzz+k),')'));
            end
        end
        options = optimset('MaxFunEvals',10000);
        Km_store1 = strcat('@(k)(',Km_store,')');
        Km_store2 = strcat('@(k)(-',Km_store,')');
        Km_store1 = strrep(Km_store1,'kA(1)',strcat('k(',num2str(length(kinit)-1),')'));
        Km_store2 = strrep(Km_store2,'kA(1)',strcat('k(',num2str(length(kinit)-1),')'));
        g1 = str2func(Km_store1);
        g2 = str2func(Km_store2);
        [~,MMparameters_final(i).KmvalueLB{j,2}] = fmincon(g1,kinit,[],[],[],[],klo,kup,[],options);
        MMparameters_final(i).KmvalueLB{j,1} = MM_symbolic(i).MM{j,1};
        MMparameters_final(i).KmvalueUB{j,1} = MM_symbolic(i).MM{j,1};
        [~,MMparameters_final(i).KmvalueUB{j,2}] = fmincon(g2,kinit,[],[],[],[],klo,kup,[],options);
        MMparameters_final(i).KmvalueUB{j,2} = -MMparameters_final(i).KmvalueUB{j,2};
        
        %MMparameters_final(i).KmvalueLB{j,2} = eval(MMparameters_final(i).Km{j,2});
        %k = MM_symbolic(i).kUB;
        %kI = MM_symbolic(i).kIUB;
        %MMparameters_final(i).KmvalueUB{j,1} = MM_symbolic(i).MM{j,1};
        %MMparameters_final(i).KmvalueUB{j,2} = eval(MMparameters_final(i).Km{j,2});
    end
    %%
    for j = 1:size(MM_symbolic(i).kcatterms,1)
        for k = 1:length(substrate)
            if k == 1
                MMparameters_final(i).kcat{j,1} = strrep(char(MM_symbolic(i).kcatterms{j,1}),strcat(substrate{k},'*'),'');
                MMparameters_final(i).kcat{j,1} = strrep(char(MMparameters_final(i).kcat{j,1}),strcat('*',substrate{k}),'');
            else
                MMparameters_final(i).kcat{j,1} = strrep(char(MMparameters_final(i).kcat{j,1}),strcat(substrate{k},'*'),'');
                MMparameters_final(i).kcat{j,1} = strrep(char(MMparameters_final(i).kcat{j,1}),strcat('*',substrate{k}),'');
            end
        end
        for k = 1:length(product)
            if k == 1
                MMparameters_final(i).kcat{j,1} = strrep(char(MMparameters_final(i).kcat{j,1}),strcat(product{k},'*'),'');
                MMparameters_final(i).kcat{j,1} = strrep(char(MMparameters_final(i).kcat{j,1}),strcat('*',product{k}),'');
            else
                MMparameters_final(i).kcat{j,1} = strrep(char(MMparameters_final(i).kcat{j,1}),strcat(product{k},'*'),'');
                MMparameters_final(i).kcat{j,1} = strrep(char(MMparameters_final(i).kcat{j,1}),strcat('*',product{k}),'');
            end  
        end
        for k = 1:length(inhibitor)
            if k == 1
                MMparameters_final(i).kcat{j,1} = strrep(char(MMparameters_final(i).kcat{j,1}),strcat(inhibitor{k},'*'),'');
                MMparameters_final(i).kcat{j,1} = strrep(char(MMparameters_final(i).kcat{j,1}),strcat('*',inhibitor{k}),'');
            else
                MMparameters_final(i).kcat{j,1} = strrep(char(MMparameters_final(i).kcat{j,1}),strcat(inhibitor{k},'*'),'');
                MMparameters_final(i).kcat{j,1} = strrep(char(MMparameters_final(i).kcat{j,1}),strcat('*',inhibitor{k}),'');
            end
        end
        for k = 1:length(activator)
            if k == 1
                MMparameters_final(i).kcat{j,1} = strrep(char(MMparameters_final(i).kcat{j,1}),strcat(activator{k},'*'),'');
                MMparameters_final(i).kcat{j,1} = strrep(char(MMparameters_final(i).kcat{j,1}),strcat('*',activator{k}),'');
            else
                MMparameters_final(i).kcat{j,1} = strrep(char(MMparameters_final(i).kcat{j,1}),strcat(activator{k},'*'),'');
                MMparameters_final(i).kcat{j,1} = strrep(char(MMparameters_final(i).kcat{j,1}),strcat('*',activator{k}),'');
            end
        end
        
        klo = [MM_symbolic(i).kLB;MM_symbolic(i).kILB];
        kup = [MM_symbolic(i).kUB;MM_symbolic(i).kIUB];
        kinit = (klo+kup)/2;
        %kI = MM_symbolic(i).kILB;
        zzz = length(MMparam(i).param);
        kcat_store = MMparameters_final(i).kcat{j,1};
        for k = 1:length(MM_symbolic(i).kIUB)
            kk = strcat('kI(',num2str(k),')');
            if strfind(Km_store,kk)
                strrep(Km_store,kk,strcat('k(',num2str(zzz+k),')'));
            end
        end
        options = optimset('MaxFunEvals',10000);
        kcat_store1 = strcat('@(k)(',kcat_store,')');
        kcat_store2 = strcat('@(k)(-',kcat_store,')');
        g1 = str2func(kcat_store1);
        g2 = str2func(kcat_store2);
        [~,MMparameters_final(i).kcatvalueUB{j,1}] = fmincon(g2,kinit,[],[],[],[],klo,kup,[],options);
        MMparameters_final(i).kcatvalueUB{j,1} = -MMparameters_final(i).kcatvalueUB{j,1};
        [~,MMparameters_final(i).kcatvalueLB{j,1}] = fmincon(g1,kinit,[],[],[],[],klo,kup,[],options);
        
               
        
        
        %MMparameters_final(i).kcatvalueLB{j,1} = eval(MMparameters_final(i).kcat{j,1});
        %k = MM_symbolic(i).kUB;
        %kI = MM_symbolic(i).kIUB;
        %MMparameters_final(i).kcatvalueUB{j,1} = eval(MMparameters_final(i).kcat{j,1});
    end
    
end

k = 1;
l = 1;
VMAX = {'f','b'};
for i = 1:size(MMparameters_final,2)
    for j = 1:size(MMparameters_final(i).KmvalueLB,1)
        Out_Km{k,1} = MMparam(i).id;
        Out_Km{k,2} = MMparameters_final(i).KmvalueLB{j,1};
        Out_Km{k,3} = MMparameters_final(i).KmvalueLB{j,2};
        Out_Km{k,4} = MMparameters_final(i).KmvalueUB{j,2};
        k = k + 1;
    end
    for j = 1:size(MMparameters_final(i).kcatvalueUB)
        Out_Vmax{l,1} = MMparam(i).id;
        Out_Vmax{l,2} = VMAX{j};
        Out_Vmax{l,3} = MMparameters_final(i).kcatvalueLB{j};
        Out_Vmax{l,4} = MMparameters_final(i).kcatvalueUB{j};
        l = l + 1;
    end
        
end

%Km apparent for comparison with experimental data

endcheck = 1;