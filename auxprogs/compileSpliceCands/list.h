#ifndef _List_H
#define _List_H

#define LINELENGTH 256

struct LineData
{
			char name[LINELENGTH];
			long int startPos;
			long int endPos;
			float averageCoverage;
			char strand;
};
typedef struct LineData *PtrToLineData;
typedef PtrToLineData Data;
struct Node;
typedef struct Node *PtrToNode;
typedef PtrToNode List;
typedef PtrToNode Position;

List        MakeEmpty(List L);
int         IsEmpty(List L);
int         IsLast(Position P, List L);
void        Delete(Data X, List L);
Position    FindPrevious(Data X, List L);
void        Insert(Data X, List L, Position P);
void        DeleteList(List L);
Position    Header(List L);
Position    First(List L);
Position    Advance(Position P);
Data 	    Retrieve(Position P);

#endif    /* _List_H */

