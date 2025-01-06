# Template pour dépots github

# Suggestion d'organisation des fichiers
- `data/` Fichiers de données (brutes et après traitement)
    -  `data-raw`
    - `data-clean`
- `R` Scripts R
     - `funs` scripts pour fonctions
     - `sims` scripts pour simulations
     - `analyses` scripts pour analyses de données
- `outputs/` Sorties des scripts R
    - `fig` Figures
    - `mods` Modèles sauvegardés
        - `sims` Modèles pour simulations
        - `analysis` Modèles pour anayse des données
    - `quarto` Rapports Quarto
- `docs` Pages html pour rapports Quarto

# Instructions pour publier vos rapports
Il existe plusieurs façon de publier vos rapports reproductibles en ligne (voir [ici](https://quarto.org/docs/publishing/github-pages.html#render-to-docs)). 
Le plus simple est d’utiliser la suite d’incantations suivante directement depuis le terminal:

- `quarto render` pour (re)compiler votre rapport
- `git add .` pour ajouter les changements produits
- `git commit -m "votre message ici"` pour documenter vos changements
- `git push` pour pousser vers votre branche principale sur github

Utilisez les instructions [ici](https://quarto.org/docs/publishing/github-pages.html) pour finaliser l’affichage de la page sur github.



# Lien pour rapport Quarto
https://metatox-morico.github.io/morico-gh-template/
