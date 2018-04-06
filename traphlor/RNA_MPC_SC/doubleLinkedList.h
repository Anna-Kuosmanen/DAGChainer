/*
 * doubleLinkedList.h
 *
 * Created on June 12th, 2014
 *	Author: aekuosma
 *
 * Customized for the purpose of merging the subpath lists.
 */

#ifndef LINKEDLIST_H
#define LINKEDLIST_H

using namespace std;
template <class T>
class doubleLinkedList
{
public:
	struct node{
		T data;
		node *prev;
        	node *next;
	};
	doubleLinkedList(){
		head=tail=NULL;
	}

	doubleLinkedList(vector<T> x, int n) {
		node *q;
		node *p=new node;   //create first node
		p->data=x[0];
		p->next=NULL;
		p->prev=NULL;
		head=p;
		for(int i=1;i<x.size();i++){
			q=p;
			p=p->next=new node;
			p->data=x.at(i);
			p->next=NULL;
			p->prev=q;
		}
		tail=p;

	}

	~doubleLinkedList(){
		while(head) {
			node* temp(head);
			head = head->next;
			delete temp;
		}
	}

	node* getHead() {
		return head;
	}

	node* getTail() {
		return tail;
	}

	void push_front(T n) {
		node* p=head;
		node* q=new node;
		q.data=n;
		p->prev = q;
		head=q;
		q->next=p;
		
	}

	void push_back(T n){
		node* p=tail;
		node* q=new node;
		q.data=n;
		p->next=q;
		q->prev=p;
		tail=q;
	}

	void delete_node(T n){
		node* p=head;
		while(p->next != NULL) {
			if(p->data == n){
				node* before=p->prev;
				node* after=p->next;
				if(before != NULL) {
					before->next=after;
					if(p==tail)
						tail=before;
				}
				if(after != NULL) {
					after->prev=before;
					if(p==head)
						head=after;
				}
				delete p;
				break;
			}
		}
	}

	void printList() {
		cout << "List: ";
		node* p=head;
		while(p != NULL) {
			cout << p->data << " ";
			p = p->next;
		}
		cout << endl;
	}

	// Merges two lists
	static void mergeLists(doubleLinkedList *prefix, doubleLinkedList *suffix, int overlap_length) {
		// Delete overlap nodes from second
		for(int i=0;i<overlap_length;i++)
				suffix->delete_node(suffix->head->data);

		// Move the pointers
		suffix->head->prev=prefix->tail;
		prefix->tail->next=suffix->head;
		prefix->tail=suffix->tail;
	}

private:
	node *head;
	node *tail;

	// No copying allowed, it'll get messy with the pointers
	doubleLinkedList& operator=( const doubleLinkedList& other );

};

#endif //LINKEDLIST_H
