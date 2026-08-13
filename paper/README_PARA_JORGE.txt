=======================================================================
PAPER THPE — PAQUETE PARA ARXIV (version inglesa y version bilingue)
Preparado el 10 de agosto de 2026
Para: Jorge Ordonez Mora (revision y verificacion)
=======================================================================

CONTENIDO DE ESTE PAQUETE

  THPE_arxiv_BILINGUAL.tex   *** VERSION RECOMENDADA PARA ENVIAR ***
                             Ingles completo primero, castellano
                             completo despues, en un solo documento.
                             16 paginas, 30 referencias compartidas.
  THPE_arxiv_BILINGUAL.pdf   Compilado, para lectura.

  THPE_arxiv_EN.tex          Version solo en ingles (10 paginas), por
                             si preferis enviar unicamente esa.
  THPE_arxiv_EN.pdf          Compilado, para lectura.

  README_PARA_JORGE.txt      Este archivo.

  FALTAN LAS FIGURAS: copiar a la misma carpeta, desde
  01_THPE\REPO_GITHUB\ o desde GitHub:
      THPE_v18_dr2_resultados.png
      THPE_v18_dr2_corner.png
  Los documentos compilan igualmente sin ellas (van envueltas en
  \IfFileExists), pero los huecos de figura quedarian vacios.

-----------------------------------------------------------------------
POR QUE UNA VERSION BILINGUE (politica verificada)
-----------------------------------------------------------------------
Desde el 11 de febrero de 2026, arXiv exige que TODOS los envios
incluyan una version inglesa completa del articulo, no solo del
resumen. Simultaneamente permite y fomenta los envios multilingues: el
mismo documento puede contener el articulo integro en dos idiomas.

Requisitos oficiales que este documento ya cumple:
  1. La version INGLESA aparece PRIMERO (requisito explicito).
  2. La version en el otro idioma va despues, completa.
  3. Se indica el idioma original: hay una "Note on languages" al
     final del abstract ingles.
  4. Formato TeX admitido.

Fuente: info.arxiv.org/help/faq/multilang.html y blog de arXiv
(anuncios de 21/11/2025 y 13/01/2026).

Ventaja: el articulo queda accesible a cualquier lector hispanohablante
sin perder alcance internacional. Coherente con un trabajo realizado en
Espana por autores independientes.

Si preferis simplicidad, la version solo inglesa es igualmente valida y
es lo mas habitual en cosmologia. La decision es vuestra.

-----------------------------------------------------------------------
COMO COMPILARLO
-----------------------------------------------------------------------
Opcion facil: Overleaf (www.overleaf.com, cuenta gratuita). Subir el
.tex y los dos .png y compilar. Trae RevTeX y babel de serie.

En local (TeX Live / MiKTeX completos):
    pdflatex THPE_arxiv_BILINGUAL.tex
    pdflatex THPE_arxiv_BILINGUAL.tex     (dos veces, por referencias)

Nota: los PDF de este paquete se compilaron con la clase 'article'
porque el entorno de preparacion no tenia RevTeX. El contenido es
identico; con RevTeX saldra a dos columnas, con aspecto de Physical
Review.

Para el castellano se recomienda anadir, si no da problemas:
    \usepackage[spanish,english]{babel}
Mejora la division silabica. No es imprescindible.

-----------------------------------------------------------------------
QUE REVISAR (checklist)
-----------------------------------------------------------------------
CRITICO — sin esto no se envia:
[ ] 1. Cotejar los 13 valores de DESI DR2 y sus 6 coeficientes de
       correlacion (en THPE_fit_v18.py, lista DESI_DR2) contra la
       tabla oficial de arXiv:2503.14738. Todo el resultado descansa
       en esa transcripcion.
[ ] 2. Verificar la paginacion de las referencias pendientes: Damour &
       Esposito-Farese 1992, y las de conjuntos causales (Sorkin 1991
       y 2000; Ahmed 2004 y 2013; Barrow 2007).
[ ] 3. Comprobar la Ec. (9) (cota de Cassini): que omega_BD =
       w F/(F')^2 y el paso a |xi| < 2.5e-3/kappa sean correctos
       dimensionalmente.
[ ] 4. Validar fisicamente la Seccion 5 (refutacion de Compton).

DE CRITERIO — decisiones de autor:
[ ] 5. Destino de gamma en la Ec. (1): mantenerlo con la exclusion
       declarada (como esta) o eliminarlo del modelo base.
[ ] 6. Enviar la version bilingue o solo la inglesa.
[ ] 7. Lectura completa del ingles por si alguna frase suena forzada.
[ ] 8. Revisar el texto de la declaracion de asistencia de IA.
[ ] 9. Revisar la traduccion castellana: es fiel al ingles, pero
       conviene comprobar si algun termino tecnico admite mejor
       equivalente en espanol.

-----------------------------------------------------------------------
COMO ENVIARLO A ARXIV
-----------------------------------------------------------------------
1. Cuenta en arxiv.org (gratuita).

2. ENDORSEMENT. Para publicar por primera vez en astro-ph.CO hace
   falta que un autor ya establecido en esa categoria os avale. No es
   revision del contenido: solo certifica que el trabajo es apropiado
   para el archivo. Se solicita desde la propia web, que genera un
   codigo para enviar a esa persona.
   Sugerencia de contactos: IFT UAM-CSIC (Madrid), ICCUB (Barcelona).
   El repositorio publico con resultado negativo y controles de
   coherencia es una buena carta de presentacion.

3. Categoria primaria sugerida: astro-ph.CO
   Secundaria: gr-qc

4. Subir el .tex y los .png. arXiv compila el LaTeX en su servidor;
   NO se sube el PDF.

5. En los metadatos, indicar que es un envio multilingue y cual es el
   idioma original (castellano).

6. Licencia recomendada: CC BY 4.0, coherente con la MIT del
   repositorio de codigo.

-----------------------------------------------------------------------
NOTA SOBRE EL ENFOQUE DEL PAPER
-----------------------------------------------------------------------
Este articulo NO presenta una teoria demostrada. Presenta:
  - cotas observacionales medidas (alpha < 0.008, beta < 0.07, gamma
    excluido);
  - la localizacion teorica de la extension dentro del espacio de
    teorias ya conocido;
  - dos advertencias metodologicas de utilidad general;
  - una prediccion falsable con registro temporal.

Esa es su fortaleza y asi debe defenderse. Las secciones VI.F y VI.G y
las conclusiones delimitan explicitamente el alcance y reconocen el
precedente de Sorkin (1990), lo que refuerza la credibilidad en lugar
de restarla.

-----------------------------------------------------------------------
MATERIAL DE APOYO
-----------------------------------------------------------------------
Expediente completo (informes de campana, fichas del programa teorico,
bibliografia verificada, cronologia): carpeta FINAL PROYECTO THPE.
Codigo e historial de correcciones: github.com/moimora1/THPE (v1.8).
=======================================================================
