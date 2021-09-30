function traceIsocritereTailleHessien(beta2, beta1, taille, age, Hessien, grad, beta2Current, beta1Current, numFig)

% Fonction qui trace les isovaleurs du critere d'erreur quadratique
% moyenne en fonction des parametres beta1 et beta2.
%
% _ beta1 = vecteur des valeurs prises par le parametre beta1 dans le modele
% _ beta2 = vecteur des valeurs prises par le parametre beta2 dans le modele
% _ taille = vecteur contenant la taille de l'enfant
% _ age = vecteur contenant l'age de l'enfant
% _ Hessien = matrice Hessienne du critere au point courant
% _ grad = vecteur du gradient du critere au point courant
% _ beta1Current = valeur du parametre beta1 a l'iteration courrante
% _ beta2Current = valeur du parametre beta2 a l'iteration courrante
% _ numFig = numero de la figure sur laquelle sont tracees les isovaleurs du
% critere


taille = taille(:);
age = age(:);

dbeta1 = beta1 - beta1Current;
dbeta2 = beta2 - beta2Current;


J = (taille - beta1Current*(1+age).^beta2Current).' * (taille - beta1Current*(1+age).^beta2Current);

isoValeur = zeros(length(dbeta1),length(dbeta2));

for ind1 = 1:length(dbeta2)
	for ind2 = 1:length(dbeta1)
		isoValeur(ind2, ind1) = J + grad.' * [dbeta1(ind2);dbeta2(ind1)] + 1/2 * [dbeta1(ind2);dbeta2(ind1)].' * Hessien * [dbeta1(ind2);dbeta2(ind1)];
	end
end

valeur = linspace(min(min(isoValeur)), J, 5);

figure(numFig);
	hold on

	[c,h] = contour(beta2, beta1, isoValeur, valeur); hold on
	set(h,'LineStyle',':','color','r','LineWidth',2)
	
		







