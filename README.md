# PPI Project
Scripts de experimentos sobre las PPIs.

Descargar PPI-Project y en la consola ejecutar lo siguiente:

```
bash deepwalk.sh
```
El mismo ejecutará lo siguiente:
* *"packages.py"* (instalará python packages que se usarán)
* *"process.py"* (ejecutará entre otras cosas Deep Walk, t-SNE y DBSCAN)
* *"analize_network.py"* (generará archivos con parámetros de las redes y distribuciones de grado)

Al finalizar la ejecución, dentro de "PPI-Project" se crearán nuevas carpetas con archivos, entre ellas "proteins" que contendrá archivos .txt por especie dentro de los cuales aparecen clusters con proteinas.

