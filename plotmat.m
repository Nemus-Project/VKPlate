close all
clear all
clc
%%
ofast=readtable("benchmark_biharm_cpp_-Ofast.csv");
gp=table2array(ofast(:,1));
Ofast=table2array(ofast(:,2));
%
ot=readtable("benchmark_biharm_cpp_1_-O3.csv");

Ot=table2array(ot(:,2));
%
matl=readtable("benchmark_biharm_matlab_1.csv");

Matl=table2array(matl(:,2));
%
pyth=readtable("benchmark_biharm_python_1.csv");

Pyth=table2array(pyth(:,2));

%
eig1=readtable("benchmark_eigs_cpp_1_-O3.csv");
gpeig1=table2array(eig1(:,1));
Eig1=table2array(eig1(:,2));
%
%
eig2=readtable("benchmark_eigs_matlab_1.csv");
gpeig2=table2array(eig2(:,1));
Eig2=table2array(eig2(:,2));

%%

ungp=unique(gp);

minof=ones(length(ungp),1)*1e236;
minot=ones(length(ungp),1)*1e236;
minmat=ones(length(ungp),1)*1e236;
minpyt=ones(length(ungp),1)*1e236;
mineigc=ones(length(ungp),1)*1e236;
mineigmat=ones(length(ungp),1)*1e236;

for i= 1:length(unique(gp))

    testgp=gp-ungp(i);

    for j=1:length(gp)

        if testgp(j)==0
            
            if testgp(j)==0
        
                if Ofast(j) < minof(i)

                    minof(i)=Ofast(j);

                end



                if Ot(j) < minot(i)

                    minot(i)=Ot(j);

                end

                if Pyth(j) < minpyt(i)

                    minpyt(i)=Pyth(j);

                end

                if Matl(j) < minmat(i)

                    minmat(i)=Matl(j);

                end

            end
        end
    end
end

%%
figure
plot(ungp,minof,"LineWidth",4)
hold on
plot(ungp,minot,"LineWidth",4)
plot(ungp,minpyt,"LineWidth",4)
plot(ungp,minmat,"LineWidth",4)

legend("Ofast","Ot","Python","Matlab")
xlabel("Gridpoints number")
ylabel("Time (s)")
set(gca,"FontSize",23)

%%

figure
plot(gpeig1,Eig1,"LineWidth",4)
hold on
plot(gpeig2,Eig2,"LineWidth",4)


legend("C","Matlab")
xlabel("Gridpoints number")
ylabel("Time (s)")
set(gca,"FontSize",23)