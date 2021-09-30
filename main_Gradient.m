% Script pour optimiser le critere par la methode de descente du gradient a
% pas fixe


clear all
close all
clc

% valeurs relevees
taille	= [0.55; 0.63; 0.66; 0.73; 0.80; 0.82; 0.86; 0.93]; % vecteur colonne des valeurs relevees de la taille
age		= [0; 2; 3; 6; 10; 12; 16; 24]; % vecteur colonne des valeurs relevees de l'age




% Parametres
rho		= 0.01;		
beta1_0 = 0.4;
beta2_0 = 0.38;
nbItMax = 20;
method = "levenberg";


% 1. Calcul des isovaleurs du critere (voir l'aide de la fonction traceIsocritereTaille)
beta1 = 0.25:1/200:0.75;
beta2 = 0:1/200:0.5;
traceIsocritereTaille(beta2, beta1, taille, age, 2)




% 2. Descente de gradient
	% a. Initialisation
		beta = [beta1_0 ; beta2_0];		% valeur initiale des parametres
		J(1) = 1/length(taille)*sum(abs((taille-beta1_0*(1+age).^beta2_0).^2));		% Valeur du critere pour les parametres initiaux
        critere = [J];
        normegrad = [sqrt(((1/length(taille))*sum((-(1+age).^beta2_0).*2.*(taille-beta1_0.*(1+age).^beta2_0)))^2+((1/length(taille))*sum((-beta1_0.*log(1+age).*(1+age).^beta2_0).*2.*(taille-beta1_0.*(1+age).^beta2_0)))^2)];
		
		% Affichage de la valeur courante des parametres sur la figure des
		% isovaleur du critere (A COMPLETER)

		
	% b. Iterations
		for ind = 2:nbItMax
			% Calcul du gradient
				gradJ = [(1/length(taille))*sum((-(1+age).^beta(2, ind-1)).*2.*(taille-beta(1, ind-1).*(1+age).^beta(2, ind-1)));
                    (1/length(taille))*sum((-beta(1, ind-1).*log(1+age).*(1+age).^beta(2, ind-1)).*2.*(taille-beta(1, ind-1).*(1+age).^beta(2, ind-1)))];
            % Calcul du hessien
           
                 H = [(1/length(taille))*sum((-(1+age).^beta(2, ind-1)).*2.*(-(1+age).^beta(2, ind-1))), (1/length(taille))*sum((-log(1+age).*(1+age).^beta(2, ind-1).*2.*(taille-beta(1, ind-1).*(1+age).^beta(2, ind-1)))+(-(1+age).^beta(2, ind-1)).*2.*(-beta(1, ind-1).*log(1+age).*(1+age).^beta(2, ind-1)));
                      (1/length(taille))*sum((-log(1+age).*(1+age).^beta(2, ind-1).*2.*(taille-beta(1, ind-1).*(1+age).^beta(2, ind-1)))+(-(1+age).^beta(2, ind-1)).*2.*(-beta(1, ind-1).*log(1+age).*(1+age).^beta(2, ind-1))), (1/length(taille))*sum((-beta(1, ind-1).*(log(1+age).^2).*(1+age).^beta(2, ind-1)).*2.*(taille-beta(1, ind-1).*(1+age).^beta(2, ind-1))+(-beta(1, ind-1).*log(1+age).*(1+age).^beta(2, ind-1)).*2.*(-beta(1, ind-1).*log(1+age).*(1+age).^beta(2, ind-1)))
                     ];
			% mise a jour des parametres
            
            if method == "pasvar"
                beta(:,ind) = beta(:,ind-1) -rho*gradJ;
                delta = 1/length(taille)*sum(abs((taille-beta(1,ind)*(1+age).^beta(2,ind)).^2)) - 1/length(taille)*sum(abs((taille-beta(1, ind-1)*(1+age).^beta(2, ind-1)).^2));
                if delta > 0
                    rho = 0.5*rho;
                    beta(:,ind) = beta(:,ind-1);
                elseif delta <= 0
                    rho = 2*rho;
                end
            elseif method == "newton"
                beta(:, ind) = beta(:,ind-1) -inv(H)*gradJ;
            elseif method == "levenberg"
                beta(:, ind) = beta(:,ind-1) -inv(H+rho*eye(2))*gradJ;
                delta = 1/length(taille)*sum(abs((taille-beta(1,ind)*(1+age).^beta(2,ind)).^2)) - 1/length(taille)*sum(abs((taille-beta(1, ind-1)*(1+age).^beta(2, ind-1)).^2));
                if delta > 0
                    rho = 0.5*rho;
                elseif delta <= 0
                    rho = 2*rho;
                end
            else
                beta(:,ind) = beta(:,ind-1) -rho*gradJ;
            end
            
            critere(ind) = 1/length(taille)*sum(abs((taille-beta(1, ind)*(1+age).^beta(2, ind)).^2));
            normegrad(ind) = sqrt(gradJ(1)^2+ gradJ(2)^2);
            
			% Affichage des courbes
                figure(1)
                hold on
                subplot(4,1,1)
                plot(1:ind, beta(1,:))
                title('Evolution de beta1 en fonction du nombre d''itération')
                subplot(4,1,2);
                plot(1:ind, beta(2,:));
                title('Evolution de beta2 en fonction du nombre d''itération')
                subplot(4,1,3);
                plot(1:ind, critere(:));
                title('Evolution du critère en fonction du nombre d''itération')
                subplot(4,1,4);
                plot(1:ind, normegrad(:));
                title('Evolution de la norme du gradient en fonction du nombre d''itération')
               
                hold off
                figure(2);
                hold on
                plot(beta(2, ind), beta(1,ind), 'r*')
                title('Evolution de beta par rapport aux isovaleurs du critère')
                hold off
                
                %traceIsocritereTaille(beta2, beta1, taille, age, 2)
               
        end
        
        figure(3)
        hold on
        
        x = 0:1/100:24;
        y = beta(1, nbItMax)*(1+x).^(beta(2, nbItMax));
   
        
        plot(age, taille, '-o');
        plot(x, y);
        title('Modèle obtenu par rapport aux points de mesure')
       
  
        
    
		
		

