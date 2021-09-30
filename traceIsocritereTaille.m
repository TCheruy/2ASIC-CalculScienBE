function traceIsocritereTaille(beta2, beta1, taille, age, numFig)

% Fonction qui trace les isovaleurs du critere d'erreur quadratique
% moyenne en fonction des parametres beta1 et beta2.
%
% _ beta1 = vecteur des valeurs prises par le parametre beta1 dans le modele pour le trace de la courbe
% _ beta2 = vecteur des valeurs prises par le parametre beta2 dans le modele pour le trace de la courbe
% _ taille = vecteur contenant la taille mesuree de l'enfant
% _ age = vecteur contenant l'age ou la taille de l'enfant a ete mesuree
% _ numFig = numero de la figure sur laquelle sont tracees les isovaleurs du
% critere

age     = age(:);
taille  = taille(:);

isoValeur = zeros(length(beta1),length(beta2));

nbPts = length(age);

for ind1 = 1:length(beta2)
	for ind2 = 1:length(beta1)
		isoValeur(ind2, ind1) = 1/nbPts * (taille - beta1(ind2)*(1+age).^beta2(ind1)).' * (taille - beta1(ind2)*(1+age).^beta2(ind1));
	end
end

valeur = exp(linspace(log(min(min(isoValeur))),log(max(max(isoValeur))) , 16));

figure(numFig);clf
	axs = axes;
		set(axs, 'FontSize', 20)

	contour(beta2, beta1, isoValeur, valeur); hold on
	
	grid on
	
	
	colorbar
		







