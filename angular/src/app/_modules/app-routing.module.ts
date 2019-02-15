import { NgModule } from '@angular/core';
import { CommonModule } from '@angular/common';
import { RouterModule, Routes } from '@angular/router';

import { IntroductionComponent, AnalysisComponent,
    HelpComponent, FaqComponent, DemoComponent } from '../_components/index';


const appRoutes: Routes = [
    // Component: 导航到此路由时，路由器需要创建的组件（HeroesComponent）
    { path: 'analysis', component: AnalysisComponent },
    { path: 'index', component: IntroductionComponent },
    { path: '', redirectTo: '/index' , pathMatch: 'full' },
    { path: 'help', component: HelpComponent },
    { path: 'faq', component: FaqComponent },
    { path: 'demo', component: DemoComponent },
  ];
  
  
  @NgModule({
    imports: [
      CommonModule,
      RouterModule.forRoot(appRoutes)
    ],
    // 把RouterModule添加到路由模块的exports中，
    // 以便关联模块（比如AppModule）中的组件可以访问路由模块中的声明，比如RouterLink 和 RouterOutlet。
    exports: [RouterModule],
    // declarations: []  // 无declarations！声明是关联模块的任务。
  })
  export class AppRoutingModule { }