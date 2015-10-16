#include <stdlib.h>
#include <stdio.h>
#include "list.h"
#include <string.h>

/* Define Linked List Node */
struct Node 
{
	Data Element;
	Position    Next;
};


List MakeEmpty(List L) 
{
	if (L != NULL) 
	{
		DeleteList(L);
  	}
	L = malloc( sizeof(struct Node) );
	if (L == NULL) 
	{
		FatalError ("MakeEmpty Error: Unable to allocate more memory.");
	}
	L->Next = NULL;
	L->Element = NULL;
	return L;
}

int FatalError (char* ErrorMessage) 
{
	printf ("%s\n", ErrorMessage);
	exit;
}

/* Return "true" if L is empty */
int IsEmpty(List L) 
{
	return L->Next == NULL;
}

/* Return "true" if P is the last position in list L.
(Parameter L is unused in this implementation) */
int IsLast(Position P, List L) 
{
	return P->Next == NULL;
}



/* Delete a node from the passed in list L.
 Cell pointed to by P->Next is deleted.  */
void Delete(Data X, List L) 
{
	Position P, TmpCell;
	P = FindPrevious(X, L);
	/* Assumption use of a header node */
	if (!IsLast(P, L)) 
	{
		/* Position (X) is found; delete it */
		TmpCell = P->Next;    
		/* bypass deleted cell */
		P->Next = TmpCell->Next;
		free(TmpCell->Element);
   	 	free(TmpCell);
	}
}

/* Find the previous item of X (maybe use a double linked list instead?) */
Position FindPrevious(Data X, List L) 
{
	Position P;
	P = L;
	while (P->Next != NULL && P->Next->Element != X) 
	{
		P = P->Next;
	}
	return P;
}


/* Insert after position P */
void Insert(Data X, List L, Position P ) 
{
	Position TmpCell;
	TmpCell = malloc( sizeof( struct Node ) );
	if (TmpCell == NULL) {
    		FatalError("Out of space!");
	}
	TmpCell->Element=X;
	TmpCell->Next = P->Next;
	P->Next = TmpCell;
}


/* Delete and free complete list */
void DeleteList(List L) 
{
	Position P, Tmp;
	P = L->Next;
	L->Next = NULL;
	while(P != NULL) 
	{
		Tmp = P->Next;
		free(P->Element);
		free(P);
		P = Tmp;
	}
	free(L);
}

Position Header(List L) 
{
	return L;
}

Position First(List L) 
{
	return L->Next;
}

Position Advance(Position P) 
{
	return P->Next;
}

Data Retrieve(Position P) 
{
	return P->Element;
}


