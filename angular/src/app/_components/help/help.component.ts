// created by Xin Xiong<onlybelter@outlook.com>, https://github.com/OnlyBelter

import { Component, OnInit } from '@angular/core';
// import {ActivatedRoute} from '@angular/router';

@Component({
  selector: 'app-help',
  templateUrl: './help.component.html',
  styleUrls: ['./help.component.css']
})

export class HelpComponent implements OnInit {
  private fragment: string;

  constructor() { }

  ngOnInit() {
    // this.route.fragment.subscribe(fragment => { this.fragment = fragment; });
  }

  // ngAfterViewInit(): void {
  //   try {
  //     document.querySelector('#' + this.fragment).scrollIntoView();
  //   } catch (e) { }
  // }


}
