/*****************************************************************************
 * Houghtool - The Software Package for efficiency measuring and visualization
 * of the HT and it's variants for line detection
 *
 * Lappeenranta University of Technology, Department of Information Technology
 * Laboratory of Information Processing
 *
 * Authors: Petri Hirvonen, Jouni Ikonen (Jouni.Ikonen@lut.fi)
 *	    Pekka Kultanen, Heikki K{lvi{inen (Heikki.Kalviainen@lut.fi)
 *
 * File:    ht_visilog.h
 * Purpose: header include for read/write images having (VISILOG) format
 * Date:    Jun 1, 1993
 *****************************************************************************/

/* COM */
/* **************************************************************************
                Constantes fichier ou structures
************************************************************************** */

#define LNAMEMAX 20
        /* longueur maximum du nom + 1 */
#define LFILEMAX 80
        /* longueur maximum du nom file + 1 */
/* parametres de longueur du descripteur */

#define DESC0   256     /* en-tete de fichiers */
#define DESCIM (sizeof (struct desc_im))
#define HDRDESC (DESC0 - DESCIM)
#define PRMPHY 11

/* **************************************************************************
                Definition de la structure IMAGE 
************************************************************************** */

struct image {
/* Parametres physiques de l'image */
        /* le format de l'enveloppe physique */
        long gx;        /* Nb de points par ligne */
        long gy;        /* Nb de lignes par plan */
        long gz;        /* Nb de plans */
        long maille;    /* maillage, 0 pour rect, 1 pour hexa */        
        /* les coordonnees de l'origine de reference */
        long ox;
        long oy;        /* Coordonnees en X, Y et Z */
        long oz;
        /* Le format arithmetique */
        long gcode;     /* Flag de type de codage */
        long bstoc;     /* Nb de bits de stockage & packing*/
        long ocode;     /* Nb d octets par point apres depacking */
        long usr_hdr;
/* Les parametres associes au fichier et/ou core */
        long fd;                /* file descriptor si >0, incore sinon */
        long f_typ;
        unsigned char *ibuf;    /* adresse du premier pixel courant incore */
        long n_row;     /* numero ligne courante */
        long byte_pnt;  /* ptr sur byte (file) ou pixel (incore) courant */
        long d_offset;  /* offset du 1er octet de donnees sur data file */
        long hors_tout; /* Nombre d'octets par lignes physiques d'image */
        char name[LNAMEMAX];
        char file[LFILEMAX];
        };

/* **************************************************************************
                structure descripteur systeme fichier
************************************************************************** */

#define I_MAGICD 0x00006931
        /* magic number "\0\0i1"  en-tete  et data */

struct desc_im {
        long i_magic;   /* magic number */
/* Taille physique de l'enveloppe */
        long i_dimx;    /* nb de pixels/ligne */
        long i_dimy;    /* nb de lignes par plan */
        long i_dimz;    /* nb de plans */
        long i_res0;    /* nb d'elements par pixel */
        long i_res00;   /* flag de type de rangement des donnees */
        long i_maille;  /* type de maillage (rectang. ou hexag. ) */
        long i_res1;    /* reserve */
/* Codage arithmetique sur fichier */
        long i_code;    /* type de codage (voir image.gcode) */
        long i_bstoc;   /* nb de bits de stockage de l'element */
        long i_res2;    /* reserve */
        long i_orix;    /* Coordonnee en x de l'origine */
        long i_oriy;    /* Coordonnee en y de l'origine */
        long i_oriz;    /* Coordonnee en z de l'origine */
        long i_res3;    /* reserve */
/* Parametres de decodage de l'entete */
        long hdr_sys;   /* taille totale en octets du descripteur principal */
        long hdr_usr;   /* reserve */
        long i_res4;    /* reserve */
        long hdr_data;  /* taille en octets de l'en-tete du fichier data */
        } ;

/* **************************************************************************
                        Flags de codage des donnees
************************************************************************** */

/* Les formats standards, dans l'ordre:
        flottant
        entier 32 bits signe
        entier 16 bits signe
        entier 8 bits signe
        entier 8 bits non signe
        binaire non packe
        binaire packe */
#define I_FLOAT 0
#define I_LONG 0x01
#define I_SHORT 0x02
#define I_CHR 0x04
#define I_CHRU 0x14
#define I_LABEL 0x42
#define I_BIN 0x22
#define I_BINP 0x20
/* les formats exotiques, packes, non packes */
#define I_PACKED 0x80

/* **************************************************************************
                Constantes de divers flags
************************************************************************** */

/* Rangement, maillage, etc... */
#define I_VECT  0
#define I_RECTA 0
#define I_HEXAG 1

/* Flags de type d'acces au fichier donnees */
#define FL_STDIO 4
#define FL_PIPE 0x10
#define FL_BD 0x20
#define FL_MEM 0x40
#define FL_TEMP 8
#define FL_READ 1
#define FL_WR 2
#define FL_RDWR 3

#define R_XPND 1
#define R_BKW 0x80
/* Protections des fichiers image */
#define PMODE   0644
