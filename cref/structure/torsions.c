/*************************************************************************

   Program:    torsions
   File:       torsions.c
   
   Version:    V1.2aSTANDALONE
   Date:       13.01.97
   Function:   Generate a complete set of backbone torsion angles for a 
               protein.
   
   Copyright:  (c) Andrew C.R. Martin 1994-7
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure and Modelling Unit,
               Department of Biochemistry and Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44(0)1372 275775
               (Work) +44(0)171 419 3284
   EMail:      martin@biochem.ucl.ac.uk
               
**************************************************************************

   This code is protected by copyright.

   The routines: main(), ShowTorsions() and Usage() may be freely copied
   and distributed for no charge providing this header is included.  The
   code for these routines may be modified as required, but any
   modifications must be documented so that the person responsible can be
   identified. If someone else breaks this code, I don't want to be blamed
   for code that does not work!

   All other routines come from the library known as Bioplib and may not
   be used outside the context of this program, torsions.c. However, a
   licence for use in other programs may be obtained from the author,
   Dr. Andrew C.R.  Martin. More details of Bioplib and the licence
   agreement may be obtained from the WWW using the URL
   http://ww.biochem.ucl.ac.uk/~martin and reading the section on
   Libraries. 

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Notes:
   ======
   
**************************************************************************

   Revision History:
   =================
   V1.0  10.06.94 Original
   V1.1  16.08.94 Added -m option for Martin's output format
   V1.1STANDALONE  02.06.95 Extracted support routines from library
   V1.2STANDALONE  12.06.95 Added -c option for pseudo-CA torsions     
   V1.2aSTANDALONE 13.01.97 New versions of doReadPDB() and fsscanf() 
                            which cope better with blank lines

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>

/************************************************************************/
/* Defines
*/
#define NUMAAKNOWN 24

/* From <bioplib/MathType.h>                                            */
typedef double REAL;
#ifndef PI
#define PI (4.0 * atan(1.0))
#endif

/* From <bioplib/SysDefs.h>                                             */
#ifndef SYS_TYPES_H     /* Unix: <sys/types.h>, MS-DOS: <sys\types.h>   */
#ifndef _TYPES_         /* Ditto                                        */
typedef short           BOOL;
#endif
#endif
#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

/* From <bioplib/macros.h>                                              */
#define NEXT(x) (x)=(x)->next
#define FREELIST(y,z)   while((y)!=NULL) \
                        {  z *_freelist_macro_q; \
                           _freelist_macro_q = (y)->next; \
                           free((char *)(y)); \
                           (y) = _freelist_macro_q; \
                        }
#define ALLOCNEXT(x,y) do { (x)->next=(y *)malloc(sizeof(y));\
                         if((x)->next != NULL) { (x)->next->next=NULL; }\
                         NEXT(x); } while(0)
#define INIT(x,y) do { x=(y *)malloc(sizeof(y)); \
                    if(x != NULL) x->next = NULL; } while(0)
#define KILLLEADSPACES(y,x)                                         \
                 do \
                 {  for((y)=(x); *(y) == ' ' || *(y) == '\t'; (y)++) ; } \
                 while(0)


/* From <bioplib/pdb.h>                                                 */
#define MAXPARTIAL 8
#define SMALL      0.000001

typedef struct pdb_entry
{
   REAL x,y,z,occ,bval;
   struct pdb_entry *next;
   int  atnum;
   int  resnum;
   char junk[8];
   char atnam[8];
   char resnam[8];
   char insert[8];
   char chain[8];
}  PDB;

#define SELECT(x,w) (x) = (char *)malloc(5 * sizeof(char)); \
                    if((x) != NULL) strncpy((x),(w),5)

/************************************************************************/
/* Globals
*/
BOOL gPDBPartialOcc;
BOOL gPDBMultiNMR;
static char sTab1[]    = {'A','C','D','E','F',
                          'G','H','I','K','L',
                          'M','N','P','Q','R',
                          'S','T','V','W','Y',
                          'E','B','Z','X'
                         };
static char sTab3[][8] = {"ALA ","CYS ","ASP ","GLU ","PHE ",
                          "GLY ","HIS ","ILE ","LYS ","LEU ",
                          "MET ","ASN ","PRO ","GLN ","ARG ",
                          "SER ","THR ","VAL ","TRP ","TYR ",
                          "PCA ","ASX ","GLX ","UNK "
                         };

/************************************************************************/
/* Prototypes
*/
PDB *ReadPDB(FILE *fp, int *natom);
PDB *doReadPDB(FILE *fp, int  *natom, BOOL AllAtoms, int OccRank, 
               int ModelNum);
PDB *SelectAtomsPDB(PDB *pdbin, int nsel, char **sel, int *natom);
static BOOL StoreOccRankAtom(int OccRank, PDB multi[MAXPARTIAL], 
                             int NPartial, PDB **ppdb, PDB **pp, 
                             int *natom);
static char *FixAtomName(char *name);
REAL phi(REAL xi, REAL yi, REAL zi, REAL xj, REAL yj, REAL zj,
         REAL xk, REAL yk, REAL zk, REAL xl, REAL yl, REAL zl);
char throne(char *three);
int fsscanf(char *buffer, char *format, ...);
void CopyPDB(PDB *out, PDB *in);
int chindex(char *string, char ch);
void padterm(char *string, int  length);
int main(int argc, char **argv);
void ShowTorsions(FILE *out, PDB *pdb, REAL *tors, BOOL Radians,
                  BOOL MartinFormat);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for converting a PDB file to torsions.

   10.06.94 Original   By: ACRM
   16.08.94 Added -m option
   12.06.95 Added -c option
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   PDB  *fullpdb,
        *pdb,
        *p,
        *p1,
        *p2,
        *p3,
        *p4;
   int  natoms,
        TorNum;
   char *sel[4];
   REAL tors[3];
   BOOL Radians      = FALSE,
        MartinFormat = FALSE,
        CATorsions   = FALSE;

   argc--;
   argv++;
   
   /* Handle any switches                                               */
   if(argc)
   {
      while(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            Usage();
            return(0);
            break;
         case 'r':
            Radians = TRUE;
            break;
         case 'm':
            MartinFormat = TRUE;
            break;
         case 'c':
            CATorsions = TRUE;
            break;
         default:
            Usage();
            return(1);
            break;
         }
         argc--;
         argv++;
      }
   }

   if(argc > 2) 
   {
      Usage();
      return(1);
   }

   /* Handle any filenames                                              */
   if(argc)
   {
      if((in=fopen(argv[0],"r"))==NULL)
      {
         fprintf(stderr,"Unable to open input file: %s\n",argv[0]);
         return(1);
      }
      
      argc--;
      argv++;

      if(argc)
      {
         if((out=fopen(argv[0],"w"))==NULL)
         {
            fprintf(stderr,"Unable to open output file: %s\n",argv[0]);
            return(1);
         }
      }
   }
   
   /* Read in the structure                                             */
   if((fullpdb = ReadPDB(in, &natoms))==NULL)
   {
      fprintf(stderr,"No atoms read from PDB file\n");
      return(1);
   }

   /* Set up the atom selection and select them                         */
   SELECT(sel[0],"CA  ");
   if(!CATorsions)
   {
      SELECT(sel[1],"N   ");
      SELECT(sel[2],"C   ");
   }
   
   if((pdb = SelectAtomsPDB(fullpdb,(CATorsions?1:3),sel,&natoms))==NULL)
   {
      fprintf(stderr,"Unable to select backbone atoms from PDB \
file (no memory?)\n");
      return(1);
   }

   /* Print title                                                       */
   if(CATorsions)
   {
      fprintf(out,"Res_N    CA_N--CA_(N+1)\n");
   }
   else
   {
      
      if(MartinFormat)
         fprintf(out,"Residue    PHI      PSI      OMEGA\n");
      else
         fprintf(out,"               PHI      PSI     OMEGA\n");
   }

   fprintf(out,"--------------------------------------\n");
   

   /* Walk the linked list and calculate torsions                       */
   tors[0] = tors[1] = tors[2] = 9999.0;
   p1      = p2      = p3      = p4      = NULL;

   if(CATorsions)
   {
      for(p=pdb; p!=NULL; NEXT(p))
      {
         p1 = p2;
         p2 = p3;
         p3 = p4;
         p4 = p;
         if(p1 && p2 && p3 && p4)   /* Got all 4 atoms                  */
         {
            tors[0] = phi(p1->x, p1->y, p1->z,
                          p2->x, p2->y, p2->z,
                          p3->x, p3->y, p3->z,
                          p4->x, p4->y, p4->z);
            if(!Radians) tors[0] *= 180.0 / PI;
            fprintf(out,"   %c    %8.3f\n",throne(p2->resnam),tors[0]);
         }
         else if(p1==NULL && p2==NULL && p3==NULL && p4) /* Got 1 atom  */
         {
            fprintf(out,"   %c        -\n",throne(p4->resnam));
         }
      }
      /* Finish off by printing the last 2 residues                     */
      fprintf(out,"   %c        -\n",throne(p3->resnam));
      fprintf(out,"   %c        -\n",throne(p4->resnam));
   }
   else
   {
      for(p=pdb; p!=NULL; NEXT(p))
      {
         if(!strncmp(p->atnam,"C   ",4))
         {
            ShowTorsions(out, p, tors, Radians, MartinFormat);
            tors[0] = tors[1] = tors[2] = 9999.0;
         }
         
         /* Get pointers to four atoms in sequence                      */
         p1 = p;
         p2 = p->next;
         if(p2 != NULL) p3 = p2->next;
         if(p3 != NULL) p4 = p3->next;
         
         if(p1==NULL || p2==NULL || p3==NULL || p4==NULL)
         {
            ShowTorsions(out, p, tors, Radians, MartinFormat);
            break;
         }
         
         if(!strncmp(p->atnam,"N   ",4))
            TorNum = 1;
         else if(!strncmp(p->atnam,"CA  ",4))
            TorNum = 2;
         else if(!strncmp(p->atnam,"C   ",4))
            TorNum = 0;
         
         tors[TorNum] = phi(p1->x, p1->y, p1->z,
                            p2->x, p2->y, p2->z,
                            p3->x, p3->y, p3->z,
                            p4->x, p4->y, p4->z);
      }
   }
   return 0;
}

/************************************************************************/
/*>void ShowTorsions(FILE *out, PDB *pdb, REAL *tors, BOOL Radians,
                     BOOL MartinFormat)
   ----------------------------------------------------------------
   Input:   FILE    *out          Output file
            PDB     *pdb          PDB record pointer
            REAL    *tors         Array of torsion angles
            BOOL    Radians       Should output be in radians
            BOOL    MartinFormat  Output in Martin's required format

   Displays the torsion angles converting to degrees if Radians flag
   not set.

   10.06.94 Original    By: ACRM
   16.08.94 Added MartinFormat
*/
void ShowTorsions(FILE *out, PDB *pdb, REAL *tors, BOOL Radians,
                  BOOL MartinFormat)
{
   if(!Radians)
   {
      int i;

      for(i=0; i<3; i++)
      {
         if(tors[i] < (REAL)9990.0)
            tors[i] *= (REAL)180.0/PI;
      }
      
   }
   
   if(MartinFormat)
   {
      fprintf(out,"   %c    ", throne(pdb->resnam));

      if(tors[0] > (REAL)9998.0)
      {
         fprintf(out,"    -    %8.3f %8.3f\n",tors[1],tors[2]);
      }
      else if(tors[1] > (REAL)9998.0)
      {
         fprintf(out,"%8.3f     -        -   \n",tors[0]);
      }
      else
      {
         fprintf(out,"%8.3f %8.3f %8.3f\n",tors[0], tors[1], tors[2]);
      }
   }
   else
   {
      fprintf(out,"%5d%c %-4s %8.3f %8.3f %8.3f\n",
              pdb->resnum, pdb->insert[0], pdb->resnam,
              tors[0], tors[1], tors[2]);
   }
}
   
/************************************************************************/
/*>void Usage(void)
   ----------------
   Displays a usage message

   10.06.94 original   By: ACRM
   16.08.94 Added -m
   12.06.95 Added -c
*/
void Usage(void)
{
   fprintf(stderr,"\ntorsions V1.2. (c) 1994-5 Andrew Martin, UCL. \
Freely Distributable\n");
   fprintf(stderr,"Generates a set of backbone torsions from a PDB \
file.\n\n");
   fprintf(stderr,"Usage: torsions [-h] [-r] [<in.pdb> [<out.tor>]]\n");
   fprintf(stderr,"       -h   This help message\n");
   fprintf(stderr,"       -r   Give results in radians\n");
   fprintf(stderr,"       -m   Output format required by Martin \
Reczko\n");
   fprintf(stderr,"       -c   Generate CA-CA pseudo-torsions\n");
   fprintf(stderr,"I/O is to stdin/stdout if unspecified.\n\n");
}

/************************************************************************/
/*>PDB *ReadPDB(FILE *fp, int *natom)
   ----------------------------------
   Input:   FILE     *fp      A pointer to type FILE in which the
                              .PDB file is stored.
   Output:  int      *natom   Number of atoms read. -1 if error.
   Returns: PDB      *pdb     A pointer to the first allocated item of
                              the PDB linked list

   Reads a PDB file into a PDB linked list

   08.07.93 Written as entry for doReadPDB()
   09.07.93 Modified to return pointer to PDB
   17.03.94 Modified to handle OccRank
   06.03.95 Added value for NMR model to read (1 = first)
*/
PDB *ReadPDB(FILE *fp,
             int  *natom)
{
   return(doReadPDB(fp, natom, TRUE, 1, 1));
}

/************************************************************************/
/*>PDB *doReadPDB(FILE *fp, int *natom, BOOL AllAtoms, int OccRank,
                  int ModelNum)
   ----------------------------------------------------------------
   Input:   FILE     *fp      A pointer to type FILE in which the
                              .PDB file is stored.
            BOOL     AllAtoms TRUE:  ATOM & HETATM records
                              FALSE: ATOM records only
            int      OccRank  Occupancy ranking
            int      ModelNum NMR Model number (0 = all)
   Output:  int      *natom   Number of atoms read. -1 if error.
   Returns: PDB      *pdb     A pointer to the first allocated item of
                              the PDB linked list

   Reads a PDB file into a PDB linked list. The OccRank value indicates
   occupancy ranking to read for partial occupancy atoms.
   If any partial occupancy atoms are read the global flag 
   gPDBPartialOcc is set to TRUE.

   04.11.88 V1.0  Original
   07.02.89 V1.1  Ignore records which aren't ATOM or HETATM
   28.03.90 V1.2  Altered field widths to match PDB standard better
                  See notes above for deviations
   28.06.90 V1.2a Buffer size increased to 85 chars.
   15.02.91 V1.2b Changed comment header to match new standard.
   07.01.92 V1.3  Ignores blank lines properly
   11.05.92 V1.4  Check on EOF in while() loop, memset() buffer. 
                  ANSIed.
   01.06.92 V1.5  Documented for autodoc
   19.06.92 V1.6  Corrected use of stdlib
   01.10.92 V1.7  Changed to use fgets()
   10.06.93 V1.9  Returns 0 on failure rather than exiting
                  Replaced SIZE with sizeof(PDB) directly
   17.06.93 V2.0  Rewritten to use fsscanf()
   08.07.93 V2.1  Split from ReadPDB()
   09.07.93 V2.2  Modified to return pointer to PDB. Rewrote allocation
                  scheme.
   17.03.94 V2.3  Handles partial occupancies
                  Sets natom to -1 if there was an error to distinguish 
                  from no atoms.
                  Handles atom names which start in column 13 rather
                  than column 14. This is allowed in the standard, but
                  very rare.
                  Sets flag for partials.
   06.04.94 V2.4  Atom names starting in column 13 have their first
                  character moved to the end if it is a digit.
   03.10.94 V2.5  Check residue number as well as atom name when running
                  through alternative atoms for partial occupancy
                  Moved increment of NPartial, so only done if there
                  is space in the array. If OccRank is 0, all atoms are
                  read regardless of occupancy.
   06.03.95 V2.7  Added value for NMR model to read (0 = all)
                  No longer static. Sets gPDBMultiNMR if ENDMDL records
                  found.
   13.01.97 V2.8  Added check on return from fsscanf. Blank lines used
                  to result in duplication of the previous line since
                  fsscanf() does not reset the variables on receiving
                  a blank line. Also fixed in fsscanf().
*/
PDB *doReadPDB(FILE *fp,
               int  *natom,
               BOOL AllAtoms,
               int  OccRank,
               int  ModelNum)
{
   char     junk[8],
            atnambuff[8],
            *atnam,
            resnam[8],
            chain[4],
            insert[4],
            buffer[160],
            CurAtom[8],
            CurIns;
   int      atnum,
            resnum,
            CurRes,
            NPartial,
            ModelCount = 1;
   double   x,y,z,
            occ,
            bval;
   PDB      *pdb  = NULL,
            *p,
            multi[MAXPARTIAL];   /* Temporary storage for partial occ   */

   *natom         = 0;
   CurAtom[0]     = '\0';
   NPartial       = 0;
   gPDBPartialOcc = FALSE;
   gPDBMultiNMR   = FALSE;

   while(fgets(buffer,159,fp))
   {
      if(ModelNum != 0)          /* We are interested in model numbers  */
      {
         if(!strncmp(buffer,"ENDMDL",6))
         {
            ModelCount++;
         }

         if(ModelCount < ModelNum)   /* Haven't reached the right model */
            continue;
         else if(ModelCount > ModelNum)    /* Gone past the right model */
            break;
      }

      if(!strncmp(buffer,"ENDMDL",6))
         gPDBMultiNMR   = TRUE;
      
      if(fsscanf(buffer,"%6s%5d%1x%5s%4s%1s%4d%1s%3x%8lf%8lf%8lf%6lf%6lf",
                 junk,&atnum,atnambuff,resnam,chain,&resnum,insert,
                 &x,&y,&z,&occ,&bval) != EOF)
      {
         if((!strncmp(junk,"ATOM  ",6)) || 
            (!strncmp(junk,"HETATM",6) && AllAtoms))
         {
            /* Fix the atom name accounting for start in column 13 or 14*/
            atnam = FixAtomName(atnambuff);
            
            /* Check for full occupancy. If occupancy is 0.0 assume that 
               it is actually fully occupied; the column just hasn't been
               filled in correctly
               
               04.10.94 Read all atoms if OccRank is 0
               */
            if(occ > (double)0.999 || 
               occ < (double)SMALL || 
               OccRank == 0)
            {
               /* Trim the atom name to 4 characters                    */
               atnam[4] = '\0';
               
               if(NPartial != 0)
               {
                  if(!StoreOccRankAtom(OccRank,multi,NPartial,&pdb,&p,
                                       natom))
                  {
                     if(pdb != NULL) FREELIST(pdb, PDB);
                     *natom = (-1);
                     return(NULL);
                  }
                  
                  /* Set partial occupancy counter to 0                 */
                  NPartial = 0;
               }
               
               /* Allocate space in the linked list                     */
               if(pdb == NULL)
               {
                  INIT(pdb, PDB);
                  p = pdb;
               }
               else
               {
                  ALLOCNEXT(p, PDB);
               }
               
               /* Failed to allocate space; free up list so far & return*/
               if(p==NULL)
               {
                  if(pdb != NULL) FREELIST(pdb, PDB);
                  *natom = (-1);
                  return(NULL);
               }
               
               /* Increment the number of atoms                         */
               (*natom)++;
               
               /* Store the information read                            */
               p->atnum  = atnum;
               p->resnum = resnum;
               p->x      = (REAL)x;
               p->y      = (REAL)y;
               p->z      = (REAL)z;
               p->occ    = (REAL)occ;
               p->bval   = (REAL)bval;
               p->next   = NULL;
               strcpy(p->junk,   junk);
               strcpy(p->atnam,  atnam);
               strcpy(p->resnam, resnam);
               strcpy(p->chain,  chain);
               strcpy(p->insert, insert);
            }
            else   /* Partial occupancy                                 */
            {
               /* Set flag to say we've got a partial occupancy atom    */
               gPDBPartialOcc = TRUE;
               
               /* First in a group, store atom name                     */
               if(NPartial == 0)
               {
                  CurIns = insert[0];
                  CurRes = resnum;
                  strncpy(CurAtom,atnam,8);
               }
               
               if(strncmp(CurAtom,atnam,strlen(CurAtom)-1) || 
                  resnum != CurRes || 
                  CurIns != insert[0])
               {
                  /* Atom name has changed 
                     Select and store the OccRank highest occupancy atom
                     */
                  if(!StoreOccRankAtom(OccRank,multi,NPartial,&pdb,&p,
                                       natom))
                  {
                     if(pdb != NULL) FREELIST(pdb, PDB);
                     *natom = (-1);
                     return(NULL);
                  }
                  
                  /* Reset the partial atom counter                     */
                  NPartial = 0;
                  strncpy(CurAtom,atnam,8);
                  CurRes = resnum;
                  CurIns = insert[0];
               }
               
               if(NPartial < MAXPARTIAL)
               {
                  /* Store the partial atom data                        */
                  multi[NPartial].atnum  = atnum;
                  multi[NPartial].resnum = resnum;
                  multi[NPartial].x      = (REAL)x;
                  multi[NPartial].y      = (REAL)y;
                  multi[NPartial].z      = (REAL)z;
                  multi[NPartial].occ    = (REAL)occ;
                  multi[NPartial].bval   = (REAL)bval;
                  multi[NPartial].next   = NULL;
                  strcpy(multi[NPartial].junk,   junk);
                  strcpy(multi[NPartial].atnam,  atnam);
                  strcpy(multi[NPartial].resnam, resnam);
                  strcpy(multi[NPartial].chain,  chain);
                  strcpy(multi[NPartial].insert, insert);
                  
                  NPartial++;
               }
            }
         }
      }
   }

   if(NPartial != 0)
   {
      if(!StoreOccRankAtom(OccRank,multi,NPartial,&pdb,&p,natom))
      {
         if(pdb != NULL) FREELIST(pdb, PDB);
         *natom = (-1);
         return(NULL);
      }
   }

   /* Return pointer to start of linked list                            */
   return(pdb);
}

/************************************************************************/
/*>static BOOL StoreOccRankAtom(int OccRank, PDB multi[MAXPARTIAL], 
                                int NPartial, PDB **ppdb, PDB **pp, 
                                int *natom)
   ----------------------------------------------------------------
   Input:   int  OccRank     Occupancy ranking required (>=1)
            PDB  multi[]     Array of PDB records for alternative atom
                             positions
            int  NPartial    Number of items in multi array
   I/O:     PDB  **ppdb      Start of PDB linked list (or NULL)
            PDB  **pp        Current position in PDB linked list (or NULL)
            int  *natom      Number of atoms read
   Returns: BOOL             Memory allocation success

   Takes an array of PDB records which represent alternative atom 
   positions for an atom. Select the OccRank'th highest occupancy and
   add this one into the PDB linked list.

   To be called by doReadPDB().

   17.03.94 Original    By: ACRM
*/
static BOOL StoreOccRankAtom(int OccRank, PDB multi[MAXPARTIAL], 
                             int NPartial, PDB **ppdb, PDB **pp, 
                             int *natom)
{
   int  i,
        j,
        IMaxOcc;
   REAL MaxOcc,
        LastOcc = (REAL)0.0;
   
   if(OccRank < 1) OccRank = 1;
   
   for(i=0; i<OccRank; i++)
   {
      MaxOcc  = (REAL)0.0;
      IMaxOcc = 0;
      
      for(j=0; j<NPartial; j++)
      {
         if(multi[j].occ >= MaxOcc)
         {
            MaxOcc  = multi[j].occ;
            IMaxOcc = j;
         }
      }
      multi[IMaxOcc].occ = (REAL)0.0;

      if(MaxOcc < (REAL)SMALL) break;
      LastOcc = MaxOcc;
   }

   /* If we ran out of rankings, take the last one to be found          */
   if(MaxOcc < (REAL)SMALL)
      MaxOcc = LastOcc;

   /* Store this atom
      Allocate space in the linked list
   */
   if(*ppdb == NULL)
   {
      INIT((*ppdb), PDB);
      *pp = *ppdb;
   }
   else
   {
      ALLOCNEXT(*pp, PDB);
   }
            
   /* Failed to allocate space; error return.                           */
   if(*pp==NULL)
      return(FALSE);
               
   /* Increment the number of atoms                                     */
   (*natom)++;
               
   /* Store the information read                                        */
   (*pp)->atnum  = multi[IMaxOcc].atnum;
   (*pp)->resnum = multi[IMaxOcc].resnum;
   (*pp)->x      = multi[IMaxOcc].x;
   (*pp)->y      = multi[IMaxOcc].y;
   (*pp)->z      = multi[IMaxOcc].z;
   (*pp)->occ    = MaxOcc;
   (*pp)->bval   = multi[IMaxOcc].bval;
   (*pp)->next   = NULL;
   strcpy((*pp)->junk,   multi[IMaxOcc].junk);
   strcpy((*pp)->atnam,  multi[IMaxOcc].atnam);
   strcpy((*pp)->resnam, multi[IMaxOcc].resnam);
   strcpy((*pp)->chain,  multi[IMaxOcc].chain);
   strcpy((*pp)->insert, multi[IMaxOcc].insert);

   /* Patch the atom name to remove the alternate letter                */
   if(strlen((*pp)->atnam) > 4)
      ((*pp)->atnam)[4] = '\0';
   else
      ((*pp)->atnam)[3] = ' ';

   return(TRUE);
}

/************************************************************************/
/*>static char *FixAtomName(char *name)
   ------------------------------------
   Input:   char  *name     Atom name read from file
   Returns: char  *         Fixed atom name (pointer into name)

   Fixes an atom name by removing leading spaces, or moving a leading
   digit to the end of the string. Used by doReadPDB()

   06.04.94 Original    By: ACRM
*/
static char *FixAtomName(char *name)
{
   char *newname;
   int  len;

   /* Default behaviour, just return the input string                   */
   newname = name;

   if(name[0] == ' ')
   {
      /* Name starts in column 14, just remove leading spaces           */
      KILLLEADSPACES(newname,name);
   }
   else      /* Name starts in column 13                                */
   {
      /* If the first character is a digit, move it to the end          */
      if(isdigit(name[0]))
      {
         if((len = chindex(name,' ')) == (-1))
         {
            /* We didn't find a space in the name, so add the character
               onto the end of the string and re-terminate
            */
            len         = strlen(name);
            newname     = name+1;
            name[len]   = name[0];
            name[len+1] = '\0';
         }
         else
         {
            /* We did find a space in the name, so put the first
               character there
            */
            newname     = name+1;
            name[len]   = name[0];
         }
      }
   }
   return(newname);
}

/************************************************************************/
/*>PDB *SelectAtomsPDB(PDB *pdbin, int nsel, char **sel, int *natom)
   -----------------------------------------------------------------
   Input:   pdbin    *PDB      Input list
            nsel     int       Number of atom types to keep
            sel      **char    List of atom types to keep
   Output:  natom    *int      Number of atoms kept
   Returns:          *PDB      Output list

   Take a PDB linked list and returns a list containing only those atom 
   types specfied in the sel array.

   To set up the list of atoms to keep, define an array of pointers 
   to char:
   e.g.     char *sel[10]
   Then define the atoms in the list thus:
            SELECT(sel[0],"N   ");
            SELECT(sel[1],"CA  ");
            SELECT(sel[2],"C   ");
            SELECT(sel[3],"O   ");
   Ensure the spaces are used!!

   N.B. The routine is non-destructive; i.e. the original PDB linked 
        list is intact after the selection process

   01.03.90 Original    By: ACRM
   28.03.90 Modified to match new version of pdb.h
   24.05.90 Fixed so the variables passed in as sel[] don't 
            *have* to be 4 chars.
   17.05.93 Modified for book. Returns BOOL.
   09.07.93 Modified to return PDB pointer. Changed allocation 
            scheme. Changed back to sel[] variables *must* be 4
            chars.
*/
PDB *SelectAtomsPDB(PDB *pdbin, int nsel, char **sel, int *natom)
{
   PDB   *pdbout  = NULL,
         *p,
         *q;
   int   i;
    
   *natom = 0;
   
   /* Step through the input PDB linked list                            */
   for(p=pdbin; p!= NULL; NEXT(p))
   {
      /* Step through the selection list                                */
      for(i=0; i<nsel; i++)
      {
         /* See if there is a match                                     */
         if(!strncmp(p->atnam,sel[i],4))
         {
            /* Alloacte a new entry                                     */
            if(pdbout==NULL)
            {
               INIT(pdbout, PDB);
               q = pdbout;
            }
            else
            {
               ALLOCNEXT(q, PDB);
            }
            
            /* If failed, free anything allocated and return            */
            if(q==NULL)
            {
               if(pdbout != NULL) FREELIST(pdbout,PDB);
               *natom = 0;
               return(NULL);
            }
            
            /* Increment atom count                                     */
            (*natom)++;
            
            /* Copy the record to the output list (sets ->next to NULL) */
            CopyPDB(q, p);
            
            break;
         }
      }
   }

   /* Return pointer to start of output list                            */
   return(pdbout);
}

/************************************************************************/
/*>REAL phi(REAL xi,REAL yi,REAL zi,REAL xj,REAL yj,REAL zj,
            REAL xk,REAL yk,REAL zk,REAL xl,REAL yl,REAL zl)
   ---------------------------------------------------------
   Input:   REAL    xi,yi,zi    Input coordinates
                    xj,yj,zj
                    xk,yk,zk
                    xl,yl,zl
   Returns: REAL                The torsion angle between the 4 atoms

   Calculates the torsion angle described by 4 sets of coordinates.

   04.03.91 Original    By: ACRM
   16.06.93 Changed float to REAL
*/
REAL phi(REAL xi,
         REAL yi,
         REAL zi,
         REAL xj,
         REAL yj,
         REAL zj,
         REAL xk,
         REAL yk,
         REAL zk,
         REAL xl,
         REAL yl,
         REAL zl)
{
   REAL xij,yij,zij,
        xkj,ykj,zkj,
        xkl,ykl,zkl,
        dxi,dyi,dzi,
        gxi,gyi,gzi,
        bi,bk,ct,
        boi2,boj2,
        z1,z2,ap,s,
        bioj,bjoi;


   /* Calculate the vectors C,B,C                                       */
   xij = xi - xj;
   yij = yi - yj;
   zij = zi - zj;
   xkj = xk - xj;
   ykj = yk - yj;
   zkj = zk - zj;
   xkl = xk - xl;
   ykl = yk - yl;
   zkl = zk - zl;

   /* Calculate the normals to the two planes n1 and n2
      this is given as the cross products:
       AB x BC
      --------- = n1
      |AB x BC|

       BC x CD
      --------- = n2
      |BC x CD|
   */
   dxi = yij * zkj - zij * ykj;     /* Normal to plane 1                */
   dyi = zij * xkj - xij * zkj;
   dzi = xij * ykj - yij * xkj;
   gxi = zkj * ykl - ykj * zkl;     /* Mormal to plane 2                */
   gyi = xkj * zkl - zkj * xkl;
   gzi = ykj * xkl - xkj * ykl;

   /* Calculate the length of the two normals                           */
   bi = dxi * dxi + dyi * dyi + dzi * dzi;
   bk = gxi * gxi + gyi * gyi + gzi * gzi;
   ct = dxi * gxi + dyi * gyi + dzi * gzi;

   boi2 = 1./bi;
   boj2 = 1./bk;
   bi   = (REAL)sqrt((double)bi);
   bk   = (REAL)sqrt((double)bk);

   z1   = 1./bi;
   z2   = 1./bk;
   bioj = bi * z2;
   bjoi = bk * z1;
   ct   = ct * z1 * z2;
   if (ct >  1.0)   ct = 1.0;
   if (ct < (-1.0)) ct = -1.0;
   ap   = acos(ct);

   s = xkj * (dzi * gyi - dyi * gzi)
     + ykj * (dxi * gzi - dzi * gxi)
     + zkj * (dyi * gxi - dxi * gyi);

   if (s < 0.0) ap = -ap;

   ap = (ap > 0.0) ? PI-ap : -(PI+ap);

   return(ap);
}

/************************************************************************/
/*>char throne(char *three)
   ------------------------
   Input:   char  *three    Three letter code
   Returns: char            One letter code

   Converts 3-letter code to 1-letter code.
   Handles ASX and GLX as X
   
   29.09.92 Original    By: ACRM
   11.03.94 Modified to handle ASX and GLX in the tables
*/
char throne(char *three)
{
   int j;

   if(three[2] == 'X')
      return('X');

   for(j=0;j<NUMAAKNOWN;j++)
      if(!strncmp(sTab3[j],three,3)) return(sTab1[j]);

   /* Only get here if the three letter code was not found              */
   return('X');
}

/************************************************************************/
/*>void CopyPDB(PDB *out, PDB *in)
   -------------------------------
   Input:   PDB  *in     Input PDB record pointer
   Output:  PDB  *out    Output PDB record pointer

   Copy a PDB record, except that the ->next is set to NULL;

   12.05.92 Original    By: ACRM
*/
void CopyPDB(PDB *out,
             PDB *in)
{
   strcpy(out->junk,   in->junk);
   strcpy(out->atnam,  in->atnam);
   strcpy(out->resnam, in->resnam);
   strcpy(out->insert, in->insert);
   strcpy(out->chain,  in->chain);
   out->atnum  = in->atnum;
   out->resnum = in->resnum;
   out->x      = in->x;
   out->y      = in->y;
   out->z      = in->z;
   out->occ    = in->occ;
   out->bval   = in->bval;
   out->next   = NULL;
}

/************************************************************************/
/*>int chindex(char *string, char ch)
   ----------------------------------
   Input:      char  *string        A string
                     ch             A character for which to search
   Returns:    int                  The offset of ch in string.
   
   Returns the offset of a character in a string. -1 if not found. This is
   used in a similar manner to strchr(), but gives an offset in the string
   rather than a pointer to the character.

   10.02.91 Original
   28.05.92 ANSIed
   06.10.93 Changed name to chindex() to avoid UNIX name clash
*/
int chindex(char  *string,
            char  ch)
{
   int count;
   
   for(count=0;count<strlen(string);count++)
      if(string[count] == ch) break;
      
   if(count >= strlen(string)) count = -1;
   
   return(count);
}

/************************************************************************/
/*>int fsscanf(char *buffer, char *format, ...)
   --------------------------------------------
   Input:   char  *buffer    Buffer from which to read information
            char  *format    Format string (like scanf() et al., but see
                             restrictions below)
   Output:  ...              Scanned output variables
   Returns: int              Number of values read (EOF if end of file or
                             no specifiers found in format string)

   Hard formatted version of sscanf(). Implements FORTRAN-like rigid
   column reading out of a string.
 
   The only parsing characters recognised are:
      %<n>f    A single precision floating point number of width <n>
      %<n>lf   A double precision floating point number of width <n>
      %<n>d    An integer of width <n>
      %<n>ld   A long integer of width <n>
      %<n>u    An unsigned of width <n>
      %<n>lu   An unsigned long of width <n>
      %<n>s    A string of width <n>
      %c       A character (of width 1)
      %<n>x    <n> spaces (like FORTRAN).
   With the exception of the %c parser, the column width, <n>,
   *must* be specified.

   Blank fields read as numbers are given a value of zero.


   17.06.93 Original    By: ACRM
   12.07.93 Added %u and %lu. Corrected %s and %c to blank rather than
            NULL strings if buffer runs out. Pads string if buffer ran
            out in the middle. Takes \n in buffer as end of string.
   24.11.95 `value' was a fixed 40 character buffer. Now changed to
            allocate a suitable number of characters as required.
   13.01.97 Previously if reading from a blank line the output variables
            were unmodified since an EOF return was done immediately.
            Now the immediate EOF return only happens if the input
            buffer is a NULL variable and the EOF on blank string is
            moved to the end such that all output variables are set to 
            zero or blank before the EOF return.
*/
int fsscanf(char *buffer, char *format, ...)
{
   va_list        ap;
   char           *FormStart,
                  *BuffStart,
                  *stop,
                  form[16],           /* Store a single formatting code */
                  *value = NULL,      /* Store an item                  */
                  *ptr,
                  type;
   int            i,
                  MaxValLength = 40,  /* Initial max value width        */
                  *IntPtr,
                  NArg     = 0,
                  width    = 0;
   BOOL           LongType = FALSE;
   double         *DblPtr;
   float          *FloatPtr;
   long           *LongPtr;
   unsigned       *UPtr;
   unsigned long  *ULongPtr;

   /* Return if line is blank                                           */
   if(!buffer) return(EOF);
   
   /* Allocate initial memory for storing a value                       */
   if((value=(char *)malloc((1+MaxValLength)*sizeof(char)))==NULL)
      return(0);

   /* Start the variable argument processing                            */
   va_start(ap, format);

   /* Intialise FormStart to the start of the format string and BuffStart
      to start of input buffer
   */
   FormStart = format;
   BuffStart = buffer;

   for(;;)
   {
      /* Flag for long variables                                        */
      LongType = FALSE;
      
      /* Find the start of a % group from the format string             */
      while(*FormStart && *FormStart != '%') FormStart++;
      if(!(*FormStart)) break;      /* Exit routine                     */
   
      /* Find the next occurence of a %                                 */
      stop = FormStart+1;
      while(*stop && *stop != '%') stop++;
   
      /* Copy these format characters into our working buffer           */
      for(i=0; FormStart != stop; i++)
         form[i] = *(FormStart++);
      form[i] = '\0';

      /* Find the type we're dealing with                               */
      ptr = form + i;
      while(*ptr == '\0' || *ptr == ' ' || *ptr == '\t') ptr--;
      type = toupper(*ptr);
      
      /* Set long flag if appropriate                                   */
      if((*(ptr-1) == 'l') || (*(ptr-1) == 'L'))
         LongType = TRUE;

      /* If it's not a character, read the width from the form string   */
      width = 0;
      if(type == 'C')
      {
         width = 1;
      }
      else
      {
         for(ptr = form+1; *ptr && isdigit(*ptr); ptr++)
         {
            width *= 10;
            width += (*ptr) - '0';
         }
      }
      
      /* See if our buffer is wide enough for this item. If not, make
         more space
      */
      if(width > MaxValLength)
      {
         if((value = (char *)realloc(value, (width+1) * sizeof(char)))
            ==NULL)
         {
            /* Unable to do allocation                                  */
            va_end(ap);
            return(0);
         }
         MaxValLength = width;
      }
      

      /* Extract width characters from the input buffer. If the input 
         buffer has run out, value will be a NULL string.
      */
      stop = BuffStart + width;
      for(i=0; *BuffStart && *BuffStart != '\n' && BuffStart != stop; i++)
         value[i] = *(BuffStart++);
      value[i] = '\0';
      
      /* Act on each type                                               */
      switch(type)
      {
      case 'F':      /* A double precision or float                     */
         if(LongType)
         {
            DblPtr = va_arg(ap, double *);
            if(sscanf(value,"%lf", DblPtr) == (-1))
               *DblPtr = (double)0.0;
         }
         else
         {
            FloatPtr  = va_arg(ap, float  *);
            if(sscanf(value,"%f",  FloatPtr) == (-1))
               *FloatPtr = (float)0.0;
         }

         break;
      case 'D':      /* An integer or long int                          */
         if(LongType)
         {
            LongPtr = va_arg(ap, long *);
            if(sscanf(value,"%ld", LongPtr) == (-1))
               *LongPtr = 0L;
         }
         else
         {
            IntPtr  = va_arg(ap, int  *);
            if(sscanf(value,"%d",  IntPtr) == (-1))
               *IntPtr = 0;
         }
         break;
      case 'U':      /* An unsigned or unsigned long                    */
         if(LongType)
         {
            ULongPtr = va_arg(ap, unsigned long *);
            if(sscanf(value,"%lu", ULongPtr) == (-1))
               *ULongPtr = 0L;
         }
         else
         {
            UPtr  = va_arg(ap, unsigned  *);
            if(sscanf(value,"%u",  UPtr) == (-1))
               *UPtr = 0;
         }
         break;
      case 'S':      /* A string                                        */
         ptr = va_arg(ap, char *);
         if(value[0])                  /* Input buffer not empty        */
         {
            *(value + width) = '\0';
            strncpy(ptr, value, width+1);

            /* If the input buffer ran out in this string, pad with 
               spaces and terminate.
            */
            if(strlen(ptr) < width) padterm(ptr, width);
         }
         else                          /* Input buffer empty            */
         {
            for(i=0; i<width; i++)
               *(ptr + i)  = ' ';
            *(ptr + width) = '\0';
         }
         break;
      case 'C':      /* A character (insert a space if buffer empty)    */
         *(va_arg(ap, char *)) = (value[0] ? value[0]: ' ');
         break;
      case 'X':      /* A column to skip                                */
         /* Fall through to default action                              */
      default:
         /* Do nothing                                                  */
         ;
      }
      
      /* If not a blank column, increment arg count                     */
      if(type != 'X') NArg++;

   }
   
   /* End variable argument parsing                                     */
   va_end(ap);

   /* Free the allocated buffer                                      */
   free(value);

   /* Return number of values read or EOF if it was a blank input       */
   if(buffer[0] == '\0' || buffer[0] == '\n') return(EOF);
   return(NArg);
}

/************************************************************************/
/*>void padterm(char *string, int  length)
   ---------------------------------------
   I/O:     char  *string   String to be padded with spaces
   Input:   int   length    Required size for string

   Pads a string with spaces to length characters, then terminates it.

   06.09.91 Original    By: ACRM   
*/
void padterm(char *string,
             int  length)
{
   int i;
   
   for(i=strlen(string); i<length; i++)
      string[i] = ' ';
   string[length] = '\0';
}
