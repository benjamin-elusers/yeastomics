# SERVER
library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

    fitdata = reactive({  propfit  })

    propsub = reactive({
        fitdata() %>% dplyr::filter( property  %in% input$properties )
    })

    output$propfit = renderPlot({
        print("show property fit")
        ggplot(evo_mpc.fit, aes(y = EVO.FULL, x = MPC)) +
            geom_ribbon(aes(ymin=.fitted.evo+0.01, ymax=Inf), fill='#BB0033',alpha=0.3)+
            geom_ribbon(aes(ymin=-Inf,ymax=.fitted.evo-0.01), fill='#00BB33',alpha=0.3)+
            geom_abline(slope=0,intercept = mean(evo_mpc.fit$EVO.FULL),size=1,col='gray',linetype=2)+
            geom_point(fill='lightgray',shape=19,alpha=0.3) +
            geom_segment(aes(xend = MPC, yend = .fitted.evo), col='gray',linetype='12',size=0.3,alpha=0.5) +
            # Overlay the protein sharing a particular properties
            geom_line(mapping=aes(y=.fitted.evo),size=2) +
            geom_segment(data=propsub(), aes(xend = MPC, yend = .fitted.evo,color=property),linetype=2,size=0.5,alpha=0.8) +
            geom_point(data=propsub(), aes(color=property), shape=19,size=1.4) + scale_color_simpsons() +
            annotate("text",y = 3, x=0.8, label='Over-estimated',col='#BB0033',vjust='inward',hjust='inward',size=8)+
            annotate("text",y = 0, x=0.8, label='Under-estimated',col='#00BB33',vjust='inward',hjust='inward',size=8)+
            xlab('Protein abundance (log10 mpc)') + ylab('Mean ER (full seq.)') + theme(legend.position = 'top',legend.direction = 'vertical')

    })

    output$propres = renderPlot({
        print("show property residuals")
        ggplot(data=evo_mpc.fit, aes(y=.resid.evo,x=MPC)) +
            geom_ribbon(aes(ymin=0.01, ymax=Inf), fill='#BB0033',alpha=0.3)+
            geom_ribbon(aes(ymin=-Inf,ymax=-0.01), fill='#00BB33',alpha=0.3)+
            geom_abline(slope=0,intercept = 0,size=1)+
            geom_point(col='gray30',shape=19,alpha=0.3) +
            geom_segment(aes(xend = MPC, yend = 0),col='gray',linetype='12',size=0.3,alpha=0.5) +
            geom_segment(data=propsub(),aes(xend = MPC, y=0, yend = .resid.evo,color=property),linetype=2,size=0.5,alpha=0.8) +
            geom_point(data=propsub(), aes(color=property), shape=19,size=1.4) + scale_color_simpsons()+
            annotate("text",y = 2, x=0.8, label='Over-estimated',col='#BB0033',vjust='inward',hjust='inward',size=8,alpha=0.5)+
            annotate("text",y = -1.2, x=0.8, label='Under-estimated',col='#00BB33',vjust='inward',hjust='inward',size=8)+
            xlab('Protein abundance (log10 mpc)') + ylab('Residuals (ER)') + theme(legend.position="none")

    })

    output$sumres = renderPlot({
        print("show sum of residuals")
        propsubres = .res_prop %>% dplyr::filter(property %in% input$properties )
        YMIN = min(-0.2,propsubres$.res_evo.SR.avg)%>% round(digits = 1)
        YMAX = max(0.2,propsubres$.res_evo.SR.avg) %>% round(digits = 1)
        ggplot(propsubres,
               aes(x = reorder(str_trunc(property,w=50), -.res_evo.SR.avg),
                   y = .res_evo.SR.avg, fill=property), group=source) +
            geom_hline(yintercept = 0,size=1,col='gray',linetype=2) +
            geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=0,ymax=Inf), fill='#BB0033',alpha=0.2)+
            geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=0), fill='#00BB33',alpha=0.2)+
            annotate("text",y = Inf, x=Inf, label='Over-estimated',col='#BB0033',vjust='inward',hjust='inward',size=5)+
            annotate("text",y = -Inf, x=Inf, label='Under-estimated',col='#00BB33',vjust='inward',hjust='inward',size=5)+
            geom_bar(stat='identity',position='dodge') +
            geom_text(aes(label=N),vjust='inward',size=8) +
            theme( axis.text = element_text(size=14)) +
            ylim(YMIN,YMAX) + scale_fill_simpsons() +  coord_flip()+  theme(legend.position="none",axis.text.y = element_text(angle=45))+
            ylab("Average residual evolutionary rate") +
            xlab('Property')
    })
})
