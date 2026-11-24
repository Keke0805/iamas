/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_YY_Y_TAB_H_INCLUDED
# define YY_YY_Y_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    JACOBIAN = 258,
    DOUBLE = 259,
    FUNCTION = 260,
    DEFVAR = 261,
    DEFFIX = 262,
    SETVAR = 263,
    SETFIX = 264,
    HESSIAN = 265,
    STOICMAT = 266,
    STOCHASTIC = 267,
    DECLARE = 268,
    INITVALUES = 269,
    EQUATIONS = 270,
    FAMILIES = 271,
    INIEQUAL = 272,
    EQNEQUAL = 273,
    EQNCOLON = 274,
    SPCPLUS = 275,
    SPCEQUAL = 276,
    FAMCOLON = 277,
    ATOMDECL = 278,
    CHECK = 279,
    CHECKALL = 280,
    REORDER = 281,
    MEX = 282,
    DUMMYINDEX = 283,
    EQNTAGS = 284,
    LOOKAT = 285,
    LOOKATALL = 286,
    MONITOR = 287,
    WRITE_ATM = 288,
    WRITE_SPC = 289,
    WRITE_MAT = 290,
    LANGUAGE = 291,
    INTFILE = 292,
    DRIVER = 293,
    INLINE = 294,
    ENDINLINE = 295,
    PARAMETER = 296,
    SPCSPC = 297,
    INISPC = 298,
    INIVALUE = 299,
    EQNSPC = 300,
    EQNSIGN = 301,
    EQNCOEF = 302,
    RATE = 303,
    SPCNR = 304,
    ATOMID = 305,
    LKTID = 306,
    MNIID = 307,
    INLCTX = 308,
    INCODE = 309,
    SSPID = 310,
    EQNLESS = 311,
    EQNTAG = 312,
    EQNGREATER = 313,
    TPTID = 314,
    USEID = 315,
    FLUX = 316,
    AUTOREDUCE = 317,
    PLSPC = 318,
    UPPERCASEF90 = 319,
    MINVERSION = 320
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 76 "scan.y" /* yacc.c:1909  */

  char str[1000];

#line 124 "y.tab.h" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_Y_TAB_H_INCLUDED  */
