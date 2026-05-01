#!/bin/bash

echo "================================="
echo "  Setup Proyecto Lane-Emden"
echo "================================="

echo "Creando entorno virtual..."
python3 -m venv entorno_ASTRO

echo "Activando entorno..."
source entorno_ASTRO/bin/activate

echo "Actualizando pip..."
pip install --upgrade pip

echo "Instalando dependencias..."
pip install -r requirements.txt

echo "================================="
echo "Setup completado correctamente"
echo "================================="
echo ""
echo "Para activar el entorno manualmente usa:"
echo "source entorno_ASTRO/bin/activate"