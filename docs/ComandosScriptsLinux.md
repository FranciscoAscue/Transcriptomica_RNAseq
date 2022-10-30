---
layout: default
title: Comandos y scripts en Linux
nav_order: 3
---

# Comandos y scripts en Linux

## Comandos basicos (Desde la terminal bash)

<p style="text-align : justify" FOND FACE="calibri">Cuando inciamos el SO de Linux, nos comunicamos con la computadora a travez una interfaz que esta en constante ejecución, a este le llamaremos <b>TERMINAL</b>, aunque sabemos que la mayoria de los SO en la actualidad tienen un entorno gráfico, para los fines de este curso nos centraremos en la terminal</p>

<p style="text-align : justify" FOND FACE="calibri">Para ejecutar un comando en general se sigue la siguiente regla de operación:</p>

INGLES

`command <options> <target> --(FILE, DIRECTORY, OTHERS)`

ESPAÑOL

`comando <opciones> <blanco> --(ARCHIVO, CARPETA, OTROS)`

## BASH PIPELINES

<p align="center" width="100%">
    <img width="37%" src="https://upload.wikimedia.org/wikipedia/commons/thumb/f/f6/Pipeline.svg/1200px-Pipeline.svg.png"> 
    <img width="60%" src="https://thumbs.gfycat.com/LikelyEmbarrassedAnophelesmosquito-size_restricted.gif"> 
</p>

### Conteo de secuencias en FastQ
``` bash
fgrep -i "@S" file.fq | wc -l

cat file.fq | echo $((`wc -l`/4))
```
### Nombre de los genes anotados 
``` bash
grep $'\tgene\t' sequence.gff3 | perl -ne '/ID=([^;]+)/ and printf("%s\n", $1)'

grep $'\tgene\t' sequence.gff3 | awk '{print $9}' | cut -d';' -f1 | sed "s/ID=//g"
```
### Mesclando scripts:

``` bash
grep $'\tgene\t' sequence.gff3 | awk '{print $5-$4";"$9}'| sed 's/Name=//g' | awk -F';' '{print $3"\t"$1}' | sort -k2n  > Tamaño_genes.txt

awk '{print $3}' sequence.gff3 | sort -d | uniq -c
```

## BASH SCRIPTS

***she bang***

```bash
#!interpreter [optional-arg]

#!/bin/sh
#!/bin/bash
#!/usr/bin/python3
#!/usr/bin/awk -f
```

```bash
#!/bin/bash

##CONSTANTES

VAR1="UUU"
VAR2="UUC"

##EJECUCION

for i in {A,U,G,C}{A,U,G,C}{A,U,G,C}
do
  if [ "$i" = "$VAR1" ] || [ "$i" = "$VAR2" ]; then
    echo "$i Phe"
  fi
done

```
### Permisos

``` bash 

# cambia permisos a archivos y carpetas

# d rwx rwx rwx / (propietario) (grupos) (demas)
# r 4
# w 2		
# x 1 

chmod 775 ACGT.bash

chmod +rwx ACGT.bash

# R recursivo 

chmod -R 755 

``` 

### Correr un script

```bash
<lenguaje> script.ext

### Sin permisos

bash ACGT.bash

### Con permisos

./ACGT.bash 

### Segundo plano

./ACGT.bash&

```

### Correr en segundo plano o screen

```bash

# Añadir el caracter &
bash scrip.bash&

### Activando un screen

screen -S name

./ACGT.bash

## salir del screen CTRL A + D

## ingresar al screen nuevamente

scren -ls 

screen -r name  
```

### Detener procesos

```bash

# muestrar los procesos actuales

ps 
pstree
top -u user

#detener procesos por comando
kill PID

top -u user 
# presiona "k" y indica el PID 

htop 
# buscar el proceso y presionar "k" 
```

## SYNTAXIS

```bash
#!/bin/bash

###COMENTARIOS

###VARIABLES (numericos (enteros o flotantes), strings (""), booleanos (T / F)

VAR1=""
VAR2=""
VAR3=""

###EXECUTION ( se puede ejecutar comandos que se usan en la terminal )

echo "started at ´date´" ## ´date´ (permite ejecutar comandos dentro de un print

mkdir -p ${VAR} ## puede intercambiar variables para usarlos dentro de los comandos.


### sintaxis de comparacion

== is equal to ; if [ "$a" == "$b" ]

!= is not equal to ; if [ "$a" != "$b" ]

< is less than; if [[ "$a" < "$b" ]] 

< is greater than if [[ "$a" > "$b" ]]

AND (&&) y OR(||)

-eq is equal to if [ "$a" -eq "$b" ]

-ne is no equal to if [ "$a" -ne "$b" ]

-gt is greater then if [ "$a" -gt "$b" ]

-ge is greater then or egual to if [ "$a" -ge "$b" ]

-lt is less than if [ "$a" -lt "$b" ]

-le is less than or equal to if [ "$a" -le "$b" ]

### estructura basica de una condicional if 

if [ "$i" != "UGA" ]; then
        echo "$i phe"
fi

### estructura basica de un while

while [ number -gt 4 ]; do
       echo "loop1"
       echo "loop2"
       echo "loop3"
done

### estrucutura basica de un for

for i in a b c
do
        comand <options> target
done
```





## 3. Programas de edición en Linux

<p style="text-align : justify">Dentro de Linux tenemos programas que podemos usar para acceder a archivos, donde podemos editarlo a gusto, muchos de estos editores fueron con los que se crearon la mayoria de programas que hoy conocemos, la programación de gran parte de las aplicaciones en varios lenguajes se ha hecho y aun puede hacerse con los editores de texto</p>

### 3.1 Editor VIM

<img align="center" width="100%" src="https://media.geeksforgeeks.org/wp-content/uploads/20210131172236/intro.png">

`Vim` es un editor mejorado, basado en el editor `Vi`, para el curso se usara con mas detalle este editor. 
Para usar este editor se debe primero llamar en la terminal con el comando: 

`vim file.txt (opcional)` si existe un archivo con el nombre `file.txt` entrara para editarlo sino creara uno con ese nombre

Ahora haremos una lista de los principales comandos que se usa en vim:

`i` para empezar a escribir, tambien funciona con la tecla `Insert`, para escribir aplica los atajos de teclado que conocemos

`Esc` para poder regresar a las funciones de `vim`

`dd` para borrar una linea completa

`yy` para copiar una linea completa

`p` para pegar la linea copiada

`:q` cierra el documento

`:q!` fuerza cerrar el documento sin guardar cambios

`:wq` guarda los cambios y cierra el documento

`/patron` permite la busqueda de un patron
` n ` avanzar DELANTE en los patrones encontrados
` N ` avanzar ATRAS en los patrones encontrados

`:tabedit` permite abrir otro documento al mismo timpo `gt` con este comando puede cambiar entre ventanas

#### otros opcionales

`split` fracciona el editor en 2 ventanas

`vsplit` fracciona el editor en 2 ventanas verticalmente

`Ctrl + w + (h-j-k-l)` permite saltar entre las ventanas


### 3.2 Editor NANO

El editor `nano` muestra las opciones una ves se inicia, para iniciar es similar que `vim`:

`nano file.txt (opcional)` 

Se antepone la tecla `Ctrl` en las opciones que aparecen en el editor:

### 3.3 Editor PICO

El editor `pico` es muy similar a `nano` solo falta anteponer la tecla `Ctrl` a las letras de opciones que aparecen

Para ejecutar este comando se ejecutar en la terminal :

`pico file.txt (opcional)`